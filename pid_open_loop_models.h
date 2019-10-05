/*
* @file	pid_open_loop_models.h
*
* Analytic plant models for the real measured step-response data of the open loop test.
* Models are used to calculate the estimated PID parameters.
* These parameters come from Taylor-McLaurin expansion of the model
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"

#include "levenberg_marquardt_solver.h"
#include "pid_functions.h"
#include "pid_avg_buf.h"
#include "pid_controller.h"


// should be multiple of 2
#define MAX_DATA_POINTS 700

enum PidTuneMethod {
	// 1st order
	PID_TUNE_CHR1,
	// 2nd order approximated with 1st order
	PID_TUNE_IMC2_1,
	PID_TUNE_CHR2_1,
	// 2nd order
	PID_TUNE_CHR2,
	PID_TUNE_VDG2,
	PID_TUNE_HP2,
};

// Used as an open-loop plant model for the "manual bump test" and as an input of a transfer function
class ModelOpenLoopPlant
{
public:
	ModelOpenLoopPlant() {} 

	ModelOpenLoopPlant(const double *params_) : params((double *)params_) {
	}

	double const *getParams() const {
		return params;
	}

protected:
	double *params;
};

// Based on FOPDT model approximated from Overdamped-SOPDT model
class ModelFopdtApproximatedFromSopdt : public ModelOpenLoopPlant {
public:
	ModelFopdtApproximatedFromSopdt(double *soParams) : ModelOpenLoopPlant(p) {
		double T2divT1 = soParams[PARAM_T2] / soParams[PARAM_T];
		double T2mulT1 = soParams[PARAM_T2] * soParams[PARAM_T];
		params[PARAM_K] = soParams[PARAM_K];
		params[PARAM_T] = 0.828 + 0.812 * T2divT1 + 0.172 * exp(-6.9 * T2divT1) * soParams[PARAM_T];
		params[PARAM_L] = 1.116 * T2mulT1 / (soParams[PARAM_T] + 1.208 * soParams[PARAM_T2]) + soParams[PARAM_L];
#ifdef PID_DEBUG
		printf("Params-1: K=%g T=%g L=%g\r\n", p[PARAM_K], p[PARAM_T], p[PARAM_L]);
#endif
	}

private:
	// A storage for the new 1st order params
	double p[3];
};

// Chien-Hrones-Reswick PID implementation for the 1st order model
class ModelChienHronesReswickFirstOrder : public ModelOpenLoopPlant {
public:	
	ModelChienHronesReswickFirstOrder(const double *params_) : ModelOpenLoopPlant(params_) {
	}

	float getKp() const {
		return (float)(0.6f / params[PARAM_K]);
	}
	float getKi() const {
		return (float)(1.0f / params[PARAM_T]);
	}
	float getKd() const {
		return (float)(1.0f / (0.5f * params[PARAM_L]));
	}
};

// "IMC-PID" Rivera-Morari-Zafiriou implementation for the 1st order model
// See "Panda R.C., Yu C.C., Huang H.P. PID tuning rules for SOPDT systems: Review and some new results"
class ModelImcPidFirstOrder : public ModelOpenLoopPlant {
public:
	ModelImcPidFirstOrder(const double *params_) : ModelOpenLoopPlant(params_) {
		lambda = fmax(0.25 * params[PARAM_L], 0.2 * params[PARAM_T]);
		Kc = (2 * params[PARAM_T] + params[PARAM_L]) / (2 * params[PARAM_K] * (lambda + params[PARAM_L]));
		Ti = params[PARAM_T] + 0.5 * params[PARAM_L];
		Td = params[PARAM_T] * params[PARAM_L] / (2.0 * params[PARAM_T] + params[PARAM_L]);
	}

	float getKp() const {
		return (float)Kc;
	}
	float getKi() const {
		return (float)(Kc / Ti);
	}
	float getKd() const {
		return (float)(Kc * Td);
	}

private:
	double lambda;
	double Kc, Ti, Td;
};

// Based on "IMC-Chien" model: "Chien, I.L., IMC-PID controller design - An extension."
// "Proceedings of the IFAC adaptive control of chemical processes conference, Copenhagen, Denmark, 1988, pp. 147-152."
class ModelChienHronesReswickSecondOrder : public ModelOpenLoopPlant {
public:
	ModelChienHronesReswickSecondOrder(const double *params_) : ModelOpenLoopPlant(params_) {
		Ti = params[PARAM_T] + params[PARAM_T2];
		Td = params[PARAM_T] * params[PARAM_T2] / Ti;
		lamda = fmax(0.25 * params[PARAM_L], 0.2 * Ti);
		Kc = Ti / (params[PARAM_K] * (lamda + params[PARAM_L]));
	}

	float getKp() const {
		return (float)Kc;
	}
	float getKi() const {
		return (float)(Kc / Ti);
	}
	float getKd() const {
		return (float)(Kc * Td);
	}

protected:
	double lamda;
	double Kc, Ti, Td;
};

// Basen on Van der Grinten Model (1963)
// "Step disturbance".
class ModelVanDerGrintenSecondOrder : public ModelOpenLoopPlant {
public:
	ModelVanDerGrintenSecondOrder(const double *params_) : ModelOpenLoopPlant(params_) {
		double T12 = params[PARAM_T] + params[PARAM_T2];
		Ti = T12 + 0.5 * params[PARAM_L];
		Td = (T12 * params[PARAM_L] + 2.0 * params[PARAM_T] * params[PARAM_T2]) / (params[PARAM_L] + 2.0 * T12);
		Kc = (0.5 + T12 / params[PARAM_L]) / params[PARAM_K];
	}

	float getKp() const {
		return (float)Kc;
	}
	float getKi() const {
		return (float)(Kc / Ti);
	}
	float getKd() const {
		return (float)(Kc * Td);
	}

protected:
	double lamda;
	double Kc, Ti, Td;
};

// Based on Haalman-Pemberton model: "Haalman, A.: Adjusting controllers for a deadtime process. Control Eng. 65, 71–73 (1965)"
// Suited for overdamped response and significant delay.
class ModelHaalmanPembertonSecondOrder : public ModelChienHronesReswickSecondOrder {
public:
	ModelHaalmanPembertonSecondOrder(const double *params_) : ModelChienHronesReswickSecondOrder(params_) {
		double T1divT2 = params[PARAM_T] / params[PARAM_T2];
		double LdivT2 = params[PARAM_L] / params[PARAM_T2];
		Ti = params[PARAM_T] + params[PARAM_T2];
		Kc = 2.0 * Ti / (3 * params[PARAM_K] * params[PARAM_L]);

		if (T1divT2 >= 0.1 && T1divT2 <= 1.0 && LdivT2 >= 0.2 && LdivT2 <= 1.0) {
			Td = Ti / 4.0;
		} else {
			Td = params[PARAM_T] * params[PARAM_T2] / Ti;
		}
		
	}
};

#if 0
// Based on "IMC-Rivera/Smith" model: 
// "Proceedings of the IFAC adaptive control of chemical processes conference, Copenhagen, Denmark, 1988, pp. 147-152."
class ModelRiveraSmithSecondOrder : public ModelOpenLoopPlant {
public:
	ModelRiveraSmithSecondOrder(const double *params_) : ModelOpenLoopPlant(params_) {
		Ti = params[PARAM_T] + params[PARAM_T2];
		Td = params[PARAM_T] * params[PARAM_T2] / Ti;
		lamda = fmax(0.25 * params[PARAM_L], 0.2 * Ti);
		Kc = (params[PARAM_T] + params[PARAM_T2]) / (params[PARAM_K] * (lamda + params[PARAM_L]));
	}

	float getKp() const {
		//return (float)(Kc * (1.0 + Td / Ti));
		return (float)Kc;
	}
	float getKi() const {
		return (float)(Kc / Ti);
	}
	float getKd() const {
		return (float)(Kc * Td);
	}

private:
	double lamda;
	double Kc, Ti, Td;
};
#endif

// PID auto-tune method using Chien-Hrones-Reswick rules and SOPDT model
class PidAutoTuneChrSopdt {
public:
	void addData(float v) {
		measuredData.addDataPoint(v);
	}

	// minValue, maxValue - input step values
	// stepPoint - data index where the step was done
	// maxPoint - data index where the output was saturated
	bool findPid(PidTuneMethod method, double minValue, double maxValue, double stepPoint, double maxPoint, const double *initialParams) {
		// without offset we cannot fit the analytic curve
		double offset = findOffset(minValue, maxValue, stepPoint, maxPoint);
		StepFunction stepFunc(minValue, maxValue, offset, stepPoint);

		// set initial params
		memcpy(params, initialParams, sizeof(params));

		// create and start solver
		switch (method) {
		// 1st order
		case PID_TUNE_CHR1: {
			FirstOrderPlusDelayLineFunction func(&stepFunc, measuredData.getBuf(), measuredData.getNumDataPoints());
			const int numParams1stOrder = 3;
			LevenbergMarquardtSolver<numParams1stOrder> solver(&func, params);
			iterationCount = solver.solve();
			break;
			}
		// 2nd order
		default: {
			SecondOrderPlusDelayLineOverdampedFunction func(&stepFunc, measuredData.getBuf(), measuredData.getNumDataPoints());
			const int numParams2ndOrder = 4;
			LevenbergMarquardtSolver<numParams2ndOrder> solver(&func, params);
			iterationCount = solver.solve();
			break;
			}
		}

		switch (method) {
		case PID_TUNE_CHR1: {
			ModelChienHronesReswickFirstOrder chr(params);
			pid.pFactor = chr.getKp();
			pid.iFactor = chr.getKi();
			pid.dFactor = chr.getKd();
			break;
		}
		case PID_TUNE_IMC2_1: {
			ModelFopdtApproximatedFromSopdt fo(params);
			ModelImcPidFirstOrder imc(fo.getParams());
			pid.pFactor = imc.getKp();
			pid.iFactor = imc.getKi();
			pid.dFactor = imc.getKd();
			break;
			}
		case PID_TUNE_CHR2_1: {
			ModelFopdtApproximatedFromSopdt fo(params);
			ModelChienHronesReswickFirstOrder chr(fo.getParams());
			pid.pFactor = chr.getKp();
			pid.iFactor = chr.getKi();
			pid.dFactor = chr.getKd();
			break;
			}
		case PID_TUNE_CHR2: {
			ModelChienHronesReswickSecondOrder chr(params);
			pid.pFactor = chr.getKp();
			pid.iFactor = chr.getKi();
			pid.dFactor = chr.getKd();
			break;
		}
		case PID_TUNE_VDG2: {
			ModelVanDerGrintenSecondOrder vdg(params);
			pid.pFactor = vdg.getKp();
			pid.iFactor = vdg.getKi();
			pid.dFactor = vdg.getKd();
			break;
			}
		case PID_TUNE_HP2:{
			ModelHaalmanPembertonSecondOrder vdg(params);
			pid.pFactor = vdg.getKp();
			pid.iFactor = vdg.getKi();
			pid.dFactor = vdg.getKd();
			break;
			}
		}

		pid.offset = (float)offset;

		return true;
	}

	template <int numPoints>
	static PidAccuracyMetric simulatePid(double target1, double target2, double dTime, const pid_s & pidParams, double params[4]) {
		StoredDataInputFunction<numPoints> plantInput;
		SecondOrderPlusDelayLineOverdampedFunction plant(&plantInput, nullptr, 0);
		PidParallelController pid(pidParams);
		//PidDerivativeFilterController pid(pidParams, 10);
		double target = target1;

		PidAccuracyMetric metric;
		// guess a previous state to minimize the system "shock"
		plantInput.addDataPoint((float)(target * params[PARAM_K]) + pidParams.offset);
		// simulate over time
		for (int i = 1; i < numPoints; i++) {
			// make a step in the middle
			if (i > numPoints / 2)
				target = target2;
			// "measure" the current value of the plant
			double pidInput = plant.getEstimatedValueAtPoint(i, params);
			// wait for the controller reaction
			double pidOutput = pid.getOutput(target, pidInput, dTime);
			// apply the reaction to the plant's pidInput
			plantInput.addDataPoint((float)pidOutput);
			metric.addPoint((double)i / (double)numPoints, pidInput, target);
#ifdef PID_DEBUG
			//printf("%d %g %g %g\r\n", i, pidInput, pidOutput, target);
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

	double findOffset(double minValue, double maxValue, double stepPoint, double maxPoint) {
		if (stepPoint < 0 || maxPoint <= stepPoint || maxPoint > measuredData.getNumDataPoints())
			return 0;
		// find the real 'min value' of the measured output data (before the step function goes up).
		double avgMeasuredMin = measuredData.getAveragedData(0, (int)stepPoint);
		// find the real 'max value' of the measured output data (after the output saturation).
		double avgMeasuredMax = measuredData.getAveragedData((int)maxPoint, measuredData.getNumDataPoints() - 1);

		if (avgMeasuredMax == avgMeasuredMin)
			return 0;
		// solve the system of equations and find the offset
		return (maxValue * avgMeasuredMin - minValue * avgMeasuredMax) / (avgMeasuredMax - avgMeasuredMin);
	}

protected:
	AveragingDataBuffer<MAX_DATA_POINTS> measuredData;
	double params[4];
	pid_s pid;
	int iterationCount = 0;
};

