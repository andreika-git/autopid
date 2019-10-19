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

#include "pid_functions.h"


// should be multiple of 2
#define MAX_DATA_POINTS 700

typedef enum {
	// 1st order
	PID_TUNE_CHR1,
	PID_TUNE_AUTO1,
	// 2nd order approximated with 1st order
	PID_TUNE_IMC2_1,
	PID_TUNE_CHR2_1,
	// 2nd order
	PID_TUNE_CHR2,
	PID_TUNE_VDG2,
	PID_TUNE_HP2,
	PID_TUNE_AUTO2,
} pid_tune_method_e;

// Used as an open-loop plant model for the "manual bump test" and as an input of a transfer function
class ModelOpenLoopPlant
{
public:
	ModelOpenLoopPlant() {} 

	ModelOpenLoopPlant(const double *params_) : params((double *)params_) {
	}

	double *getParams() {
		return params;
	}

	virtual float getKp() const = 0;
	virtual float getKi() const = 0;
	virtual float getKd() const = 0;

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

	// we don't really need them, because this model is just an intermediate
	virtual float getKp() const {
		return 0.0f;
	}
	virtual float getKi() const {
		return 0.0f;
	}
	virtual float getKd() const {
		return 0.0f;
	}

private:
	// A storage for the new 1st order params
	double p[3];
};

// Standard PID model: Kc * (1 + 1/(Ti*S) + Td * S)
// This class converts in into our "Parallel" form: Kp + Ki / S + Kd * S
class ModelStandard : public ModelOpenLoopPlant {
public:
	ModelStandard(const double *params_) : ModelOpenLoopPlant(params_) {
	}

	virtual float getKp() const {
		return (float)Kc;
	}
	virtual float getKi() const {
		return (float)(Kc / Ti);
	}
	virtual float getKd() const {
		return (float)(Kc * Td);
	}

protected:
	// "Standard" PID coefs
	double Kc, Ti, Td;
};

class ModelStandardIMC : public ModelStandard {
public:
	ModelStandardIMC(const double *params_) : ModelStandard(params_) {
		lambda = fmax(0.25 * params[PARAM_L], 0.2 * Ti);
	}

protected:
	// closed-loop speed of response
	double lambda;
};


// Chien-Hrones-Reswick PID implementation for the 1st order model (generic model).
class ModelChienHronesReswickFirstOrder : public ModelStandardIMC {
public:
	ModelChienHronesReswickFirstOrder(const double *params_) : ModelStandardIMC(params_) {
		double l2 = params[PARAM_L] / 2.0;
		Ti = params[PARAM_T] + l2;
		Td = params[PARAM_T] * params[PARAM_L] / (2 * params[PARAM_T] + params[PARAM_L]);
		Kc = Ti / (params[PARAM_K] * (lambda + l2));
	}
};

// Chien-Hrones-Reswick PID implementation for the 1st order model (set-point regulation).
class ModelChienHronesReswickFirstOrderSetpoint : public ModelOpenLoopPlant {
public:	
	ModelChienHronesReswickFirstOrderSetpoint(const double *params_) : ModelOpenLoopPlant(params_) {
	}

	virtual float getKp() const {
		return (float)(0.6f / params[PARAM_K]);
	}
	virtual float getKi() const {
		return (float)(1.0f / params[PARAM_T]);
	}
	virtual float getKd() const {
		return (float)(1.0f / (0.5f * params[PARAM_L]));
	}
};

// Chien-Hrones-Reswick PID implementation for the 1st order model (disturbance rejection).
class ModelChienHronesReswickFirstOrderDisturbance : public ModelOpenLoopPlant {
public:
	ModelChienHronesReswickFirstOrderDisturbance(const double *params_) : ModelOpenLoopPlant(params_) {
	}

	virtual float getKp() const {
		return (float)(0.95f / params[PARAM_K]);
	}
	virtual float getKi() const {
		return (float)(2.4f / params[PARAM_T]);
	}
	virtual float getKd() const {
		return (float)(1.0f / (0.42f * params[PARAM_L]));
	}
};

// "IMC-PID" Rivera-Morari-Zafiriou implementation for the 1st order model
// See "Panda R.C., Yu C.C., Huang H.P. PID tuning rules for SOPDT systems: Review and some new results"
class ModelRiveraMorariFirstOrder : public ModelStandardIMC {
public:
	ModelRiveraMorariFirstOrder(const double *params_) : ModelStandardIMC(params_) {
		Kc = (2 * params[PARAM_T] + params[PARAM_L]) / (2 * params[PARAM_K] * (lambda + params[PARAM_L]));
		Ti = params[PARAM_T] + 0.5 * params[PARAM_L];
		Td = params[PARAM_T] * params[PARAM_L] / (2.0 * params[PARAM_T] + params[PARAM_L]);
	}
};

// Based on "IMC-Chien" (aka Rivera/Smith) model: "Chien, I.L., IMC-PID controller design - An extension."
// "Proceedings of the IFAC adaptive control of chemical processes conference, Copenhagen, Denmark, 1988, pp. 147-152."
class ModelChienHronesReswickSecondOrder : public ModelStandardIMC {
public:
	ModelChienHronesReswickSecondOrder(const double *params_) : ModelStandardIMC(params_) {
		Ti = params[PARAM_T] + params[PARAM_T2];
		Td = params[PARAM_T] * params[PARAM_T2] / Ti;
		Kc = Ti / (params[PARAM_K] * (lambda + params[PARAM_L]));
	}
};

// Basen on Van der Grinten Model (1963)
// "Step disturbance".
class ModelVanDerGrintenSecondOrder : public ModelStandard {
public:
	ModelVanDerGrintenSecondOrder(const double *params_) : ModelStandard(params_) {
		double T12 = params[PARAM_T] + params[PARAM_T2];
		Ti = T12 + 0.5 * params[PARAM_L];
		Td = (T12 * params[PARAM_L] + 2.0 * params[PARAM_T] * params[PARAM_T2]) / (params[PARAM_L] + 2.0 * T12);
		Kc = (0.5 + T12 / params[PARAM_L]) / params[PARAM_K];
	}
};

// Based on Haalman-Pemberton model: "Haalman, A.: Adjusting controllers for a deadtime process. Control Eng. 65, 71–73 (1965)"
// Suited for overdamped response and significant delay.
class ModelHaalmanPembertonSecondOrder : public ModelStandard {
public:
	ModelHaalmanPembertonSecondOrder(const double *params_) : ModelStandard(params_) {
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

// Based on IMC-Maclaurin model:
// Lee, Y., Park, S., Lee, M., and Brosilow, C., PID controller tuning for desired closed - loop responses for SI / SO systems. AIChE J. 44, 106–115 1998.
class ModelMaclaurinSecondOrder : public ModelStandardIMC {
public:
	ModelMaclaurinSecondOrder(const double *params_) : ModelStandardIMC(params_) {
		double T1T2 = params[PARAM_T] + params[PARAM_T2];
		double L = params[PARAM_L];
		double L2 = L * L;
		double twolL = 2.0 * lambda + L;
		Ti = T1T2 - (2.0 * lambda * lambda - L2) / (2.0 * twolL);
		Kc = Ti / (params[PARAM_K] * twolL);
		Td = Ti - T1T2 + (params[PARAM_T] * params[PARAM_T2] - L2 * L / (6.0 * twolL)) / Ti;
	}
};

class ModelAutoSolver : public ModelOpenLoopPlant {
public:
	ModelAutoSolver(const ModelOpenLoopPlant *initial) : ModelOpenLoopPlant(pidParams) {
		pidParams[PARAM_Kp] = initial->getKp();
		pidParams[PARAM_Ki] = initial->getKi();
		pidParams[PARAM_Kd] = initial->getKd();
	}

	virtual float getKp() const {
		return (float)pidParams[PARAM_Kp];
	}
	virtual float getKi() const {
		return (float)pidParams[PARAM_Ki];
	}
	virtual float getKd() const {
		return (float)pidParams[PARAM_Kd];
	}

protected:
	double pidParams[numParamsForPid];
};
