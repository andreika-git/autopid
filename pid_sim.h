/*
* @file	pid_sim.h
*
* PID Controller simulator.
* Used 
*
* @date Oct 10, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"


typedef enum {
	PID_SIM_REGULATOR = 0,	// imitate load disturbance
	PID_SIM_SERVO,	// imitate setpoint change
} pid_sim_type_e;

template <int maxPoints>
class PidSimulator {
public:
	PidSimulator(pid_sim_type_e simType_, int order, double_t target1_, double_t target2_, double_t target_, double_t dTime_, double_t modelBias, const char *outputFile_) :
		simType(simType_), target1(target1_), target2(target2_), dTime(dTime_), methodOrder(order), outputFile(outputFile_),
		plantInput(1.0 / dTime_), plant1(&plantInput, nullptr, maxPoints, modelBias), plant2(&plantInput, nullptr, maxPoints, modelBias) {
		// if we don't know the target, we want to be in the middle for safety
		target = std::isnan(target_) ? ((target1 + target2) / 2.0) : target_;
		
		if (methodOrder == 1)
			plant = (LMSFunction<4> *)&plant1;
		else
			plant = &plant2;
	}

	void setModelParams(const double_t *modelParams) {
		this->modelParams = modelParams;
		//plant->calculateAllPoints(modelParams);
	}

	PidAccuracyMetric simulate(int numPoints, const pid_s & pidParams) {
		//PidParallelController pid(pidParams);
		PidIndustrialController pid(pidParams);

		plantInput.reset();
		// guess a previous state to minimize the system "shock"
		plantInput.addDataPoint((float_t)(getSetpoint(0) / modelParams[PARAM_K]) - pidParams.offset);
		// "calm down" the PID controller to avoid huge spikes at the beginning
		pid.getOutput(getSetpoint(0), plant->calcEstimatedValuesAtPoint(0, modelParams), dTime);

		stepPoint = maxPoints / 2;

		metric.reset();

		// simulate over time
		for (int i = 0; i < numPoints; i++) {
			// make a step in the middle
			double_t target = getSetpoint(i);
			// "measure" the current value of the plant (we cannot use getEstimatedValueAtPoint() because the function input is changing)
			double_t pidInput = plant->calcEstimatedValuesAtPoint(i, modelParams) + getLoadDisturbance(i);
			// wait for the controller reaction
			double_t pidOutput = pid.getOutput(target, pidInput, dTime);

			// apply the reaction to the plant's pidInput
			plantInput.addDataPoint((float_t)pidOutput);
			// don't take into account any start-up transients, we're interested only in our step response!
			if (i >= stepPoint)
				metric.addPoint((double_t)i / (double_t)numPoints, pidInput, target);
#ifdef PID_DEBUG
			if (outputFile != nullptr)
				output_csv(outputFile, (double_t)i, pidOutput, target, pidInput);
#endif
		}

		return metric;
	}

	double_t getSetpoint(int i) const {
		switch (simType) {
		case PID_SIM_SERVO:
			return (i > stepPoint) ? target2 : target1;
		default:
			// if we don't know the target, we want to be in the middle for safety
			return target;
		}
	}

	double_t getLoadDisturbance(int i) const {
		static const double_t disturb = 0.10;
		static const double_t ampl = 0.25, period = 0.37;
		double_t d = 0;
		switch (simType) {
		case PID_SIM_REGULATOR:
			// add or subtract 10% to imitate the "load"
			d += target1 * ((i > stepPoint) ? -disturb : disturb);
			// add periodic noise
			d += sin(2.0 * 3.14159265 * (double_t)i * dTime / period) * ampl;
			return d;
		default:
			return 0.0;
		}
	}

protected:
	int methodOrder;
	int stepPoint;
	const char *outputFile;
	StoredDataInputFunction<maxPoints> plantInput;
	
	LMSFunction<4> *plant;
	FirstOrderPlusDelayLineFunction plant1;
	SecondOrderPlusDelayLineOverdampedFunction plant2;
	
	PidAccuracyMetric metric;
	double_t target, target1, target2, dTime;
	pid_sim_type_e simType;
	const double_t *modelParams;
};


// A working copy of simulator and PID controller params. Used by PidAutoTune::solveModel().
// We don't want to create a simulator instance each time we call getEstimatedValueAtPoint() from PidCoefsFinderFunction.
// So we create an instance of this class and store temporary allocated data.
template <int numPoints>
class PidSimulatorFactory {
public:
	PidSimulatorFactory(pid_sim_type_e simType, int methodOrder_, double_t target1_, double_t target2_, double_t target_, double_t dTime, double_t modelBias, const pid_s & pid_) :
		sim(simType, methodOrder_, target1_, target2_, target_, dTime, modelBias, nullptr), pid(pid_) {
	}

public:
	PidSimulator<numPoints> sim;
	pid_s pid;
};


// A special function used for PID controller analytic simulation
// This method doesn't use any "magic" formulas, instead it uses Levenberg-Marquardt solver as a "brute-force".
template <int numPoints>
class PidCoefsFinderFunction : public LMSFunction<numParamsForPid> {
public:
	PidCoefsFinderFunction(PidSimulatorFactory<numPoints> *simFactory_, const double_t *modelParams_) :
		simFactory(simFactory_), modelParams(modelParams_) {
		// try different value because this function is too complicated and not smooth
		delta = 1.0e-5;
		simFactory->sim.setModelParams(modelParams_);
	}

	virtual void justifyParams(double_t *params) const {
		// todo: limit PID coefs somehow?
		//params[PARAM_Ki] = fmax(params[PARAM_Ki], 0.00001);
	}

	// Get the total number of data points
	virtual int getNumPoints() const {
		return numPoints;
	}

	// Calculate the sum of the squares of the residuals
	virtual double_t calcMerit(double_t *params) const {
		return simulate(numPoints - 1, params).getMerit();
	}

	virtual double_t getResidual(int i, const double_t *params) const {
		return simFactory->sim.getSetpoint(i) - getEstimatedValueAtPoint(i, params);
	}

	virtual double_t calcEstimatedValuesAtPoint(int i, const double_t *params) const {
		return simulate(i, params).getLastValue();
	}

	PidAccuracyMetric simulate(int i, const double_t *params) const {
		// update params
		simFactory->pid.pFactor = (float_t)params[PARAM_Kp];
		simFactory->pid.iFactor = (float_t)params[PARAM_Ki];
		simFactory->pid.dFactor = (float_t)params[PARAM_Kd];
		// simulate PID controller response
		return simFactory->sim.simulate(i + 1, simFactory->pid);
	}

protected:
	const double_t *modelParams;
	PidSimulatorFactory<numPoints> *simFactory;
};

