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
#include "pid_sim.h"
#include "output_csv.h"


/// These settings should be obtained from the measured data
class PidAutoTuneSettings {
public:
	double_t minValue, maxValue;
	double_t stepPoint, maxPoint;
	double_t timeScale;
	double_t targetValue;
};

// PID auto-tune method using different tuning rules and FOPDT/SOPDT models
class PidAutoTune {
public:
	PidAutoTune() {
		measuredData.init();
	}

	void addData(float_t v) {
		measuredData.addDataPoint(v);
	}

	// minValue, maxValue - input step values
	// stepPoint - data index where the step was done
	// maxPoint - data index where the output was saturated
	// this is not thread-safe (because of internal static vars) but anyway it works...
	bool findPid(pid_sim_type_e simType, pid_tune_method_e method, const PidAutoTuneSettings & settings, const double_t *initialParams) {
		// save current settings & simType
		this->settings = settings;
		this->simType = simType;
		// without bias we cannot fit the analytic curve
		modelBias = findModelBias();
		
		// todo: where do we get them?
		pid.minValue = 0;
		pid.maxValue = 100;

#ifdef PID_DEBUG
		printf("* modelBias=%Lg avgMin=%Lg avgMax=%Lg\r\n", (long double)modelBias, (long double)avgMeasuredMin, (long double)avgMeasuredMax);
#endif
		StepFunction stepFunc(settings.minValue, settings.maxValue, settings.stepPoint, settings.timeScale);

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
			printf("* Params0: K=%Lg T1=%Lg T2=%Lg L=%Lg\r\n", (long double)params[PARAM_K], (long double)params[PARAM_T], (long double)params[PARAM_T2], (long double)params[PARAM_L]);
#endif
		}
		else {
			// use provided initial params
			memcpy(params, initialParams, sizeof(params));
		}

		double_t merit0, merit;
#ifdef PID_DEBUG
		printf("* Solving...\r\n");
#endif
		// create and start solver
		if (methodOrder == 1) {
			// 1st order
			FirstOrderPlusDelayLineFunction func(&stepFunc, measuredData.getBuf(), measuredData.getNumDataPoints(), modelBias);
			func.justifyParams(params);
			func.calculateAllPoints(params);
			const int numParams1stOrder = 3;
			outputFunc<numParams1stOrder>("pid_func01.csv", func, stepFunc, params);
			LevenbergMarquardtSolver<numParams1stOrder> solver(&func, params);
			merit0 = func.calcMerit(params);
			iterationCount = solver.solve(minParamT);
			outputFunc<numParams1stOrder>("pid_func1.csv", func, stepFunc, params);
			merit = func.calcMerit(params);
		}
		else {
			// 2nd order or approximated 2nd->1st order
			SecondOrderPlusDelayLineOverdampedFunction func(&stepFunc, measuredData.getBuf(), measuredData.getNumDataPoints(), modelBias);
			func.justifyParams(params);
			func.calculateAllPoints(params);
			const int numParams2ndOrder = 4;
			outputFunc<numParams2ndOrder>("pid_func02.csv", func, stepFunc, params);
			LevenbergMarquardtSolver<numParams2ndOrder> solver(&func, params);
			merit0 = func.calcMerit(params);
			iterationCount = solver.solve(minParamT);
			outputFunc<numParams2ndOrder>("pid_func2.csv", func, stepFunc, params);
			merit = func.calcMerit(params);
		}

#ifdef PID_DEBUG
		printSolverResult(iterationCount, merit0, merit);
#endif

		ModelOpenLoopPlant *model = nullptr;
		switch (method) {
		case PID_TUNE_CHR1:
		case PID_TUNE_AUTO1: {
			static ModelChienHronesReswickFirstOrder chr(params);
			model = &chr;
			break;
		}
		case PID_TUNE_IMC2_1: {
			ModelFopdtApproximatedFromSopdt fo(params);
			static ModelRiveraMorariFirstOrder imc(fo.getParams());
			model = &imc;
			break;
		}
		case PID_TUNE_CHR2_1: {
			ModelFopdtApproximatedFromSopdt fo(params);
			static ModelChienHronesReswickFirstOrder chr(fo.getParams());
			model = &chr;
			break;
		}
		case PID_TUNE_CHR2:
		case PID_TUNE_AUTO2: {
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
		// round offset and period due to the firmware limitations...
		pid.offset = round(getPidOffset(model));
		pid.periodMs = round((float_t)(1000.0 / settings.timeScale));
		pid.antiwindupFreq = 0.0;
		pid.derivativeFilterLoss = 0.0;

		pid0 = pid;

		// for "automatic" methods, we try to make the coeffs even better!
		if (method == PID_TUNE_AUTO1 || method == PID_TUNE_AUTO2) {
			ModelAutoSolver autoSolver(model);
			solveModel(method, autoSolver);
			pid.pFactor = autoSolver.getKp();
			pid.iFactor = autoSolver.getKi();
			pid.dFactor = autoSolver.getKd();
		}

		return true;
	}

	int getMethodOrder(pid_tune_method_e method) {
		switch (method) {
			// 1st order
		case PID_TUNE_CHR1:
		case PID_TUNE_AUTO1:
			return 1;
			// 2nd order
		default:
			return 2;
		}
	}

	double_t const *getParams() const {
		return params;
	}

	pid_s const &getPid() const {
		return pid;
	}

	pid_s const &getPid0() const {
		return pid0;
	}

	double_t getAvgMeasuredMin() const {
		return avgMeasuredMin;
	}

	double_t getAvgMeasuredMax() const {
		return avgMeasuredMax;
	}

	double_t getModelBias() const {
		return modelBias;
	}

	// The model output is typically shifted
	double_t findModelBias() {
		if (settings.stepPoint < 0 || settings.maxPoint <= settings.stepPoint || settings.maxPoint > measuredData.getNumDataPoints())
			return 0;
		// find the real 'min value' of the measured output data (before the step function goes up).
		avgMeasuredMin = measuredData.getAveragedData(0, (int)settings.stepPoint);
		// find the real 'max value' of the measured output data (after the output saturation).
		avgMeasuredMax = measuredData.getAveragedData((int)settings.maxPoint, measuredData.getNumDataPoints() - 1);

		if (avgMeasuredMax == avgMeasuredMin)
			return 0;
		// solve the system of equations and find the bias
		return (settings.maxValue * avgMeasuredMin - settings.minValue * avgMeasuredMax) / (settings.maxValue - settings.minValue);
	}

	// If the target value is known, we can estimate the PID offset value based on the model gain and bias.
	float_t getPidOffset(ModelOpenLoopPlant *model) const {
		if (std::isnan(settings.targetValue))
			return 0;
		return (float_t)((settings.targetValue - modelBias) / model->getParams()[PARAM_K]);
	}

	// See: Rangaiah G.P., Krishnaswamy P.R. Estimating Second-Order plus Dead Time Model Parameters, 1994.
	// Also see: "Practical PID Control", p. 169
	bool findFirstOrderInitialParams2Points(double_t *params) const {
		int i0 = (int)settings.stepPoint;
		int i1 = (int)settings.maxPoint;

		double_t dy = avgMeasuredMax - avgMeasuredMin;

		double_t t[2];
		static const double_t tCoefs[] = { 0.353, 0.853 };
		for (int i = 0; i < 2; i++) {
			t[i] = getTimeDelta(measuredData.findDataAt((float_t)(avgMeasuredMin + dy * tCoefs[i]), i0, i1));
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
	bool findSecondOrderInitialParams3Points(double_t *params) const {
		int i0 = (int)settings.stepPoint;
		int i1 = (int)settings.maxPoint;

		double_t dy = avgMeasuredMax - avgMeasuredMin;

		double_t t[3];
		static const double_t tCoefs[] = { 0.14, 0.55, 0.91 };
		for (int i = 0; i < 3; i++) {
			t[i] = getTimeDelta(measuredData.findDataAt((float_t)(avgMeasuredMin + dy * tCoefs[i]), i0, i1));
			if (t[i] == 0.0)
				return false;
		}
		double_t alpha = (t[2] - t[1]) / (t[1] - t[0]);
#if 0
		// check if usable range?
		if (alpha < 1.2323 || alpha > 2.4850) {
			return false;
		}
#endif
		double_t beta = log(alpha / (2.485 - alpha));
		double_t xi = 0.50906 + 0.51743 * beta - 0.076284 * pow(beta, 2) + 0.041363 * pow(beta, 3)
			- 0.0049224 * pow(beta, 4) + 0.00021234 * pow(beta, 5);
		double_t Tcoef = 0.85818 - 0.62907 * xi + 1.2897 * pow(xi, 2) - 0.36859 * pow(xi, 3) + 0.038891 * pow(xi, 4);
		double_t Lcoef = 1.39200 - 0.52536 * xi + 1.2991 * pow(xi, 2) - 0.36859 * pow(xi, 3) + 0.037605 * pow(xi, 4);

		double_t T = Tcoef / (t[1] - t[0]);
		// we've got T and xi, and we have to solve quadratic equation to get T1 and T2:
		// T = T1*T2
		// Xi = (T1 + T2) / (T1*T2)
		double_t det = (T * xi) * (T * xi) - 4.0 * T;

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
	bool findSecondOrderInitialParamsHarriott(double_t *params) const {
		static const HarriotFunction hfunc;

		int i0 = (int)settings.stepPoint;
		int i1 = (int)settings.maxPoint;

		double_t dy = avgMeasuredMax - avgMeasuredMin;

		double_t A1 = -measuredData.getArea(i0, i1, (float_t)avgMeasuredMax) / settings.timeScale;

		double_t t73 = getTimeDelta(measuredData.findDataAt((float_t)(avgMeasuredMin + dy * 0.73), i0, i1));
		double_t tm = i0 + 0.5 * (t73 / 1.3) * settings.timeScale;
		double_t ym = measuredData.getValue((float_t)tm);
		// normalize
		double_t ymn = (ym - avgMeasuredMin) / dy;
		// sanity check?
		if (ymn < HarriotFunction::minX || ymn > HarriotFunction::maxX)
			return false;
		double_t r = hfunc.getValue(ymn);

		params[PARAM_K] = dy / (settings.maxValue - settings.minValue);
		params[PARAM_L] = A1 - t73 / 1.3;
		params[PARAM_T] = r * t73 / 1.3;
		params[PARAM_T2] = (1.0 - r) * t73 / 1.3;

		return true;
	}

	double_t getTimeDelta(int i1) const {
		if (i1 < 0)
			return -1.0;
		return (double_t)(i1 - settings.stepPoint) / settings.timeScale;
	}

	// Use automatic LM-solver to find the best PID coefs which satisfy the minimal PID metric.
	// The initial PID coefs are already calculated using the well-known CHR method (1st or 2nd order).
	bool solveModel(pid_tune_method_e method, ModelAutoSolver & model) {
		double_t merit0, merit;
		
#ifdef PID_DEBUG
		printf("* Solving for better coefs:\r\n");
#endif
		// todo: is it correct?
		double_t dTime = pid.periodMs / 1000.0;
		const int numSimPoints = 1024;
		PidSimulatorFactory<numSimPoints> simFactory(simType, getMethodOrder(method), getAvgMeasuredMin(), getAvgMeasuredMax(), settings.targetValue, dTime, modelBias, pid);
		PidCoefsFinderFunction<numSimPoints> func(&simFactory, params);
		func.justifyParams(model.getParams());
		merit0 = func.calcMerit(model.getParams());

		// now hopefully we'll find even better coefs!
		LevenbergMarquardtSolver<numParamsForPid> solver((LMSFunction<numParamsForPid> *)&func, model.getParams());
		double lambdaForPid = 10.0;
		double minDeltaForPid = 1.e-7;
		int iterationCount = solver.solve(lambdaForPid, minDeltaForPid);
		
		merit = func.calcMerit(model.getParams());
#ifdef PID_DEBUG
		printSolverResult(iterationCount, merit0, merit);
#endif
		return true;
	}
	
	template <int numParams>
	void outputFunc(const char *fname, const AbstractDelayLineFunction<numParams> & func, const StepFunction & stepFunc, double_t *params) {
#ifdef PID_DEBUG
		//func.calculateAllPoints(params);
		for (int i = 0; i < func.getNumPoints(); i++) {
			double_t v = func.getEstimatedValueAtPoint(i, params);
			double_t sv = stepFunc.getValue((float_t)i, 0);
			output_csv(fname, (double_t)i, func.getDataPoint(i), v, sv);
		}
#endif
	}

	void printSolverResult(int iterationCount, double_t merit0, double_t merit) {
		if (iterationCount > 0)
			printf("* The solver finished in %d iterations! (Merit: %Lg -> %Lg)\r\n", iterationCount, (long double)merit0, (long double)merit);
		else
			printf("* The solver aborted after %d iterations! (Merit: %Lg -> %Lg)\r\n", -iterationCount, (long double)merit0, (long double)merit);
	}


protected:
	AveragingDataBuffer<MAX_DATA_POINTS> measuredData;
	PidAutoTuneSettings settings;
	pid_sim_type_e simType;
	double_t params[4] = { 0 };
	pid_s pid;
	pid_s pid0;	// not-optimized
	int iterationCount = 0;
	double_t avgMeasuredMin, avgMeasuredMax;
	double_t modelBias;
};

