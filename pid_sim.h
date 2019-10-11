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



template <int maxPoints>
class PidSimulator {
public:
	PidSimulator(int order, double target1_, double target2_, double dTime_, bool showOutput_) : 
		target1(target1_), target2(target2_), dTime(dTime_), methodOrder(order), showOutput(showOutput_),
		plantInput(1.0 / dTime_), plant1(&plantInput, nullptr, 0, minParamT0), plant2(&plantInput, nullptr, 0, minParamT0) {
	}

	PidAccuracyMetric simulate(int numPoints, const pid_s & pidParams, const double *params) {
		LMSFunction<4> *plant;
		if (methodOrder == 1)
			plant = (LMSFunction<4> *)&plant1;
		else
			plant = &plant2;

		PidParallelController pid(pidParams);
		//PidDerivativeFilterController pid(pidParams, 10);

		plantInput.reset();
		// guess a previous state to minimize the system "shock"
		plantInput.addDataPoint((float)(getTarget(0) / params[PARAM_K]) - pidParams.offset);
		// "calm down" the PID controller to avoid huge spikes at the beginning
		pid.getOutput(getTarget(0), plant->getEstimatedValueAtPoint(0, params), dTime);

		metric.reset();

		// simulate over time
		for (int i = 0; i < numPoints; i++) {
			// make a step in the middle
			double target = getTarget(i);
			// "measure" the current value of the plant
			double pidInput = plant->getEstimatedValueAtPoint(i, params);
			// wait for the controller reaction
			double pidOutput = pid.getOutput(target, pidInput, dTime);
			// apply the reaction to the plant's pidInput
			plantInput.addDataPoint((float)pidOutput);
			metric.addPoint((double)i / (double)numPoints, pidInput, target);
#ifdef PID_DEBUG
			if (showOutput)
				output_csv((methodOrder == 1) ? "pid_test.csv" : "pid_test2.csv", (double)i, pidInput, pidOutput, target);
#endif
		}

		return metric;
	}

	double getTarget(int i) const {
		return (i > maxPoints / 2) ? target2 : target1;
	}

protected:
	int methodOrder;
	bool showOutput;
	StoredDataInputFunction<maxPoints> plantInput;
	FirstOrderPlusDelayLineFunction plant1;
	SecondOrderPlusDelayLineOverdampedFunction plant2;
	PidAccuracyMetric metric;
	double target1, target2, dTime;
};


// A working copy of simulator and PID controller params. Used by PidAutoTune::solveModel().
// We don't want to create a simulator instance each time we call getEstimatedValueAtPoint() from PidCoefsFinderFunction.
// So we create an instance of this class and store temporary allocated data.
template <int numPoints>
class PidSimulatorFactory {
public:
	PidSimulatorFactory(int methodOrder_, double target1_, double target2_, double dTime, const pid_s & pid_) :
		sim(methodOrder_, target1_, target2_, dTime, false), pid(pid_) {
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
	PidCoefsFinderFunction(PidSimulatorFactory<numPoints> *simFactory_, const double *modelParams_) :
		simFactory(simFactory_), modelParams(modelParams_) {
	}

	virtual void justifyParams(double *params) const {
		// todo: limit PID coefs somehow?
	}

	// Get the total number of data points
	virtual int getNumPoints() const {
		return numPoints;
	}

	// Calculate the sum of the squares of the residuals
	virtual double calcMerit(double *params) const {
		return simulate(numPoints - 1, params).getItae();
	}

	virtual double getResidual(int i, const double *params) const {
		return simFactory->sim.getTarget(i) - getEstimatedValueAtPoint(i, params);
	}

	virtual double getEstimatedValueAtPoint(int i, const double *params) const {
		return simulate(i, params).getLastValue();
	}

	PidAccuracyMetric simulate(int i, const double *params) const {
		// update params
		simFactory->pid.pFactor = (float)params[PARAM_Kp];
		simFactory->pid.iFactor = (float)params[PARAM_Ki];
		simFactory->pid.dFactor = (float)params[PARAM_Kd];
		// simulate PID controller response
		return simFactory->sim.simulate(i + 1, simFactory->pid, modelParams);
	}

protected:
	const double *modelParams;
	PidSimulatorFactory<numPoints> *simFactory;
};

