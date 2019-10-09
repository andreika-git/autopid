/*
* @file	pid_auto.h
*
* PID Auto tune.
*
* @date Oct 08, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"
#include <fstream>
#include <stdarg.h>

#include "levenberg_marquardt_solver.h"
#include "pid_functions.h"
#include "pid_avg_buf.h"
#include "pid_controller.h"
#include "pid_open_loop_models.h"
#include "output_csv.h"


/// These settings should be obtained from the measured data
class PidAutoTuneSettings {
public:
	double minValue, maxValue;
	double stepPoint, maxPoint;
	double timeScale;
};

// PID auto-tune method using Chien-Hrones-Reswick rules and SOPDT model
class PidAutoTuneChrSopdt {
public:
	PidAutoTuneChrSopdt() {
		measuredData.init();
	}

	void addData(float v) {
		measuredData.addDataPoint(v);
	}

	// minValue, maxValue - input step values
	// stepPoint - data index where the step was done
	// maxPoint - data index where the output was saturated
	// this is not thread-safe (because of internal static vars) but anyway it works...
	bool findPid(PidTuneMethod method, const PidAutoTuneSettings & settings, const double *initialParams) {
		// save current settings
		this->settings = settings;
		// without offset we cannot fit the analytic curve
		double offset = findOffset();
#ifdef PID_DEBUG
		printf("* offset=%g avgMin=%g avgMax=%g\r\n", offset, avgMeasuredMin, avgMeasuredMax);
#endif
		StepFunction stepFunc(settings.minValue, settings.maxValue, offset, settings.stepPoint, settings.timeScale);

		int methodOrder = getMethodOrder(method);

		// set initial params
		if (initialParams == nullptr) {
			// guess initial params
			if (methodOrder == 1) {
				findFirstOrderInitialParams2Points(params);
			} else {
				if (!findSecondOrderInitialParamsHarriott(params)) {
					if (!findSecondOrderInitialParams3Points(params)) {
						// this is a bad scenario, but we don't have a choice...
						// first, find FOPDT params (at least it doesn't fail)
						findFirstOrderInitialParams2Points(params);
						// then imitate FOPDT with SODPT, when T1+T2 = T, T1*T2 -> 0
						params[PARAM_T2] = 0.01;
						// and hope that our solver will do the rest...
					}
				}
			}
#ifdef PID_DEBUG
			printf("* Params0: K=%g T1=%g T2=%g L=%g\r\n", params[PARAM_K], params[PARAM_T], params[PARAM_T2], params[PARAM_L]);
#endif
		}
		else {
			// use provided initial params
			memcpy(params, initialParams, sizeof(params));
		}

		// try to solve several times unless we see some progress (iterationCount > 1)
		for (double minParamT = minParamT0; minParamT <= minParamT1; minParamT *= 10.0) {
			double merit0, merit;
#ifdef PID_DEBUG
			printf("* Solving with minParamT=%g...\r\n", minParamT);
#endif
			// create and start solver
			if (methodOrder == 1) {
				// 1st order
				FirstOrderPlusDelayLineFunction func(&stepFunc, measuredData.getBuf(), measuredData.getNumDataPoints(), minParamT);
				func.justifyParams(params);
				const int numParams1stOrder = 3;
				outputFunc<numParams1stOrder>("pid_func0.csv", func, stepFunc, params);
				LevenbergMarquardtSolver<numParams1stOrder> solver(&func, params);
				merit0 = solver.calcMerit(params);
				iterationCount = solver.solve(minParamT);
				outputFunc<numParams1stOrder>("pid_func.csv", func, stepFunc, params);
				merit = solver.calcMerit(params);
			}
			else {
				// 2nd order or approximated 2nd->1st order
				SecondOrderPlusDelayLineOverdampedFunction func(&stepFunc, measuredData.getBuf(), measuredData.getNumDataPoints(), minParamT);
				func.justifyParams(params);
				const int numParams2ndOrder = 4;
				outputFunc<numParams2ndOrder>("pid_func0.csv", func, stepFunc, params);
				LevenbergMarquardtSolver<numParams2ndOrder> solver(&func, params);
				merit0 = solver.calcMerit(params);
				iterationCount = solver.solve(minParamT);
				outputFunc<numParams2ndOrder>("pid_func.csv", func, stepFunc, params);
				merit = solver.calcMerit(params);
			}

#ifdef PID_DEBUG
			if (iterationCount > 0)
				printf("* The solver finished in %d iterations! Merit (%g -> %g)\r\n", iterationCount, merit0, merit);
			else
				printf("* The solver aborted after %d iterations! Merit (%g -> %g)\r\n", -iterationCount, merit0, merit);
#endif
			if (abs(iterationCount) > 1)
				break;
		}

		ModelOpenLoopPlant *model = nullptr;
		switch (method) {
		case PID_TUNE_CHR1: {
			static ModelChienHronesReswickFirstOrder chr(params);
			model = &chr;
			break;
		}
		case PID_TUNE_IMC2_1: {
			ModelFopdtApproximatedFromSopdt fo(params);
			static ModelImcPidFirstOrder imc(fo.getParams());
			model = &imc;
			break;
		}
		case PID_TUNE_CHR2_1: {
			ModelFopdtApproximatedFromSopdt fo(params);
			static ModelChienHronesReswickFirstOrder chr(fo.getParams());
			model = &chr;
			break;
		}
		case PID_TUNE_CHR2: {
			static ModelChienHronesReswickSecondOrder chr(params);
			model = &chr;
			break;
		}
		case PID_TUNE_VDG2: {
			static ModelVanDerGrintenSecondOrder vdg(params);
			model = &vdg;
			break;
		}
		case PID_TUNE_HP2: {
			static ModelHaalmanPembertonSecondOrder hp(params);
			model = &hp;
			break;
		}
		default:
			return false;
		}
		pid.pFactor = model->getKp();
		pid.iFactor = model->getKi();
		pid.dFactor = model->getKd();

		pid.offset = (float)offset;
		// todo: where do we get them?
		pid.minValue = 0;
		pid.maxValue = 100;

		return true;
	}

	int getMethodOrder(PidTuneMethod method) {
		switch (method) {
			// 1st order
		case PID_TUNE_CHR1:
			return 1;
			// 2nd order
		default:
			return 2;
		}
	}

	template <int numPoints>
	static PidAccuracyMetric simulatePid(int order, double target1, double target2, double dTime, const pid_s & pidParams, const double *params) {
		StoredDataInputFunction<numPoints> plantInput(1.0 / dTime);
		FirstOrderPlusDelayLineFunction plant1(&plantInput, nullptr, 0, minParamT0);
		SecondOrderPlusDelayLineOverdampedFunction plant2(&plantInput, nullptr, 0, minParamT0);
		LMSFunction<4> *plant;
		if (order == 1)
			plant = (LMSFunction<4> *)&plant1;
		else
			plant = &plant2;

		PidParallelController pid(pidParams);
		//PidDerivativeFilterController pid(pidParams, 10);
		double target = target1;

		PidAccuracyMetric metric;
		// guess a previous state to minimize the system "shock"
		plantInput.addDataPoint((float)(target / params[PARAM_K]) - pidParams.offset);
		// simulate over time
		for (int i = 1; i < numPoints; i++) {
			// make a step in the middle
			if (i > numPoints / 2)
				target = target2;
			// "measure" the current value of the plant
			double pidInput = plant->getEstimatedValueAtPoint(i, params);
			// wait for the controller reaction
			double pidOutput = pid.getOutput(target, pidInput, dTime);
			// apply the reaction to the plant's pidInput
			plantInput.addDataPoint((float)pidOutput);
			metric.addPoint((double)i / (double)numPoints, pidInput, target);
#ifdef PID_DEBUG
			output_csv((order == 1) ? "pid_test.csv" : "pid_test2.csv", (double)i, pidInput, pidOutput, target);
#endif
		}

		return metric;
	}


	double const *getParams() const {
		return params;
	}

	pid_s const &getPid() const {
		return pid;
	}

	double getAvgMeasuredMin() const {
		return avgMeasuredMin;
	}

	double getAvgMeasuredMax() const {
		return avgMeasuredMax;
	}

	double findOffset() {
		if (settings.stepPoint < 0 || settings.maxPoint <= settings.stepPoint || settings.maxPoint > measuredData.getNumDataPoints())
			return 0;
		// find the real 'min value' of the measured output data (before the step function goes up).
		avgMeasuredMin = measuredData.getAveragedData(0, (int)settings.stepPoint);
		// find the real 'max value' of the measured output data (after the output saturation).
		avgMeasuredMax = measuredData.getAveragedData((int)settings.maxPoint, measuredData.getNumDataPoints() - 1);

		if (avgMeasuredMax == avgMeasuredMin)
			return 0;
		// solve the system of equations and find the offset
		return (settings.maxValue * avgMeasuredMin - settings.minValue * avgMeasuredMax) / (avgMeasuredMax - avgMeasuredMin);
	}

	// See: Rangaiah G.P., Krishnaswamy P.R. Estimating Second-Order plus Dead Time Model Parameters, 1994.
	// Also see: "Practical PID Control", p. 169
	bool findFirstOrderInitialParams2Points(double *params) const {
		int i0 = (int)settings.stepPoint;
		int i1 = (int)settings.maxPoint;

		double dy = avgMeasuredMax - avgMeasuredMin;

		double t[2];
		static const double tCoefs[] = { 0.353, 0.853 };
		for (int i = 0; i < 2; i++) {
			t[i] = getTimeDelta(measuredData.findDataAt((float)(avgMeasuredMin + dy * tCoefs[i]), i0, i1));
			if (t[i] < 0.0)
				return false;
		}

		params[PARAM_K] = dy / (settings.maxValue - settings.minValue);
		params[PARAM_T] = 0.67 * (t[1] - t[0]);
		params[PARAM_L] = 1.3 * t[0] - 0.29 * t[1];

		return true;
	}

	// See: Rangaiah G.P., Krishnaswamy P.R. Estimating Second-Order plus Dead Time Model Parameters, 1994.
	// Also see: "Practical PID Control", p. 187
	bool findSecondOrderInitialParams3Points(double *params) const {
		int i0 = (int)settings.stepPoint;
		int i1 = (int)settings.maxPoint;

		double dy = avgMeasuredMax - avgMeasuredMin;

		double t[3];
		static const double tCoefs[] = { 0.14, 0.55, 0.91 };
		for (int i = 0; i < 3; i++) {
			t[i] = getTimeDelta(measuredData.findDataAt((float)(avgMeasuredMin + dy * tCoefs[i]), i0, i1));
			if (t[i] == 0.0)
				return false;
		}
		double alpha = (t[2] - t[1]) / (t[1] - t[0]);
#if 0
		// check if usable range?
		if (alpha < 1.2323 || alpha > 2.4850) {
			return false;
		}
#endif
		double beta = log(alpha / (2.485 - alpha));
		double xi = 0.50906 + 0.51743 * beta - 0.076284 * pow(beta, 2) + 0.041363 * pow(beta, 3)
			- 0.0049224 * pow(beta, 4) + 0.00021234 * pow(beta, 5);
		double Tcoef = 0.85818 - 0.62907 * xi + 1.2897 * pow(xi, 2) - 0.36859 * pow(xi, 3) + 0.038891 * pow(xi, 4);
		double Lcoef = 1.39200 - 0.52536 * xi + 1.2991 * pow(xi, 2) - 0.36859 * pow(xi, 3) + 0.037605 * pow(xi, 4);

		double T = Tcoef / (t[1] - t[0]);
		// we've got T and xi, and we have to solve quadratic equation to get T1 and T2:
		// T = T1*T2
		// Xi = (T1 + T2) / (T1*T2)
		double det = (T * xi) * (T * xi) - 4.0 * T;

		params[PARAM_K] = dy / (settings.maxValue - settings.minValue);
		if (det < 0) {
			// that's a hack :(
			// we ignore xi and equalize T and T2 to match the higher-order coefficient.
			if (T < 0)
				return false;
			if (xi > 0) {
				params[PARAM_T2] = T * xi;
				params[PARAM_T] = T / params[PARAM_T2];
			} else {
				params[PARAM_T2] = params[PARAM_T] = sqrt(T);
			}
		} else {
			params[PARAM_T2] = T * xi + sqrt(det) / 2.0;	// we take larger root
			params[PARAM_T] = T / params[PARAM_T2];
		}
		params[PARAM_L] = t[1] - T * Lcoef;

		return true;
	}

	// See: Harriott P. Process control (1964). McGraw-Hill. USA.
	// Also see: "Practical PID Control", p. 182
	bool findSecondOrderInitialParamsHarriott(double *params) const {
		static const HarriotFunction hfunc;

		int i0 = (int)settings.stepPoint;
		int i1 = (int)settings.maxPoint;

		double dy = avgMeasuredMax - avgMeasuredMin;

		double A1 = -measuredData.getArea(i0, i1, (float)avgMeasuredMax) / settings.timeScale;

		double t73 = getTimeDelta(measuredData.findDataAt((float)(avgMeasuredMin + dy * 0.73), i0, i1));
		double tm = i0 + 0.5 * (t73 / 1.3) * settings.timeScale;
		double ym = measuredData.getValue((float)tm);
		// normalize
		double ymn = (ym - avgMeasuredMin) / dy;
		// sanity check?
		if (ymn < HarriotFunction::minX || ymn > HarriotFunction::maxX)
			return false;
		double r = hfunc.getValue(ymn);

		params[PARAM_K] = dy / (settings.maxValue - settings.minValue);
		params[PARAM_L] = A1 - t73 / 1.3;
		params[PARAM_T] = r * t73 / 1.3;
		params[PARAM_T2] = (1.0 - r) * t73 / 1.3;

		return true;
	}

	double getTimeDelta(int i1) const {
		if (i1 < 0)
			return -1.0;
		return (double)(i1 - settings.stepPoint) / settings.timeScale;
	}

	
	template <int numParams>
	void outputFunc(const char *fname, const AbstractDelayLineFunction<numParams> & func, const StepFunction & stepFunc, double *params) {
#ifdef PID_DEBUG
		for (int i = 0; i < func.getNumPoints(); i++) {
			double v = func.getEstimatedValueAtPoint(i, params);
			double sv = stepFunc.getValue((float)i, 0);
			output_csv(fname, (double)i, func.getDataPoint(i), v, sv);
		}
#endif
	}

protected:
	AveragingDataBuffer<MAX_DATA_POINTS> measuredData;
	PidAutoTuneSettings settings;
	double params[4] = { 0 };
	pid_s pid;
	int iterationCount = 0;
	double avgMeasuredMin, avgMeasuredMax;
};

