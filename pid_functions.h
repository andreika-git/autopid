/*
* @file	pid_functions.h
*
* Functions used by PID tuner
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"

#include "levenberg_marquardt_solver.h"
#include "pid_avg_buf.h"

// This helps to differentiate with respect to 'delay' axis (while finding the Hessian matrix)
#define INTERPOLATED_STEP_FUNCTION

// Generic Step model params:
enum {
	PARAM_K = 0,	// K = Gain
	PARAM_T,		// T = Time constant
	PARAM_L,		// L = Delay (dead time)

	// 2nd order params
	PARAM_T2,		// T2 = Time constant

};

static const double minParamValue = 0.0001;

class InputFunction {
public:
	virtual double getValue(double i) const = 0;
};

// Heaviside step function interpolated between 'min' and 'max' values with 'stepPoint' time offset
class StepFunction : public InputFunction
{
public:
	StepFunction(double minValue_, double maxValue_, double offset_, double stepPoint_) :
				minValue(minValue_), maxValue(maxValue_), offset(offset_), stepPoint(stepPoint_) {
	}

	virtual double getValue(double i) const {
#ifdef INTERPOLATED_STEP_FUNCTION
		// the delay parameter L may not be integer, so we have to interpolate between the closest input values (near and far in the past)
		int I = (int)i;
		double fract = i - I;	// 0 = choose near value, 1 = choose far value
		// find two closest input values for the given delay
		double vNear = (I < stepPoint) ? minValue : maxValue;
		double vFar = (I + 1 < stepPoint) ? minValue : maxValue;
		// interpolate
		return offset + vFar * fract + vNear * (1.0f - fract);
#else
		return offset + ((i < stepPoint) ? minValue : maxValue);
#endif
	}

private:
	double minValue, maxValue;
	// stepPoint is float because we have AveragingDataBuffer, and the time axis may be scaled
	double stepPoint;
	// needed to use PARAM_K coefficient properly; also offset is used by PID
	double offset;
};

template <int numPoints>
class StoredDataInputFunction : public InputFunction {
public:
	void addDataPoint(float v) {
		// todo: support data scaling
		assert(inputData.getNumDataPoints() <= numPoints);
		inputData.addDataPoint(v);
	}

	virtual double getValue(double i) const {
		return inputData.getValue((float)i);
	}

private:
	AveragingDataBuffer<numPoints> inputData;
};

// Abstract indirect transfer function used for step response analytic simulation
template <int numParams>
class AbstractDelayLineFunction : public LMSFunction<numParams> {
public:
	AbstractDelayLineFunction(const InputFunction *input, const float *measuredOutput, int numDataPoints) {
		dataPoints = measuredOutput;
		inputFunc = input;
		numPoints = numDataPoints;
	}

	virtual double getResidual(int i, const double *params) const {
		return dataPoints[i] - getEstimatedValueAtPoint(i, params);
	}

	virtual double getEstimatedValueAtPoint(int i, const double *params) const = 0;

	// Get the total number of data points
	virtual double getNumPoints() {
		return numPoints;
	}

protected:
	const InputFunction *inputFunc;
	const float *dataPoints;
	int numPoints;
};

// FODPT indirect transfer function used for step response analytic simulation.
// Used mostly as an approximate model for chemical processes?
// The Laplace representation is: K * exp(-L*s) / (T*s + 1)
class FirstOrderPlusDelayLineFunction : public AbstractDelayLineFunction<3> {
public:
	FirstOrderPlusDelayLineFunction(const InputFunction *input, const float *measuredOutput, int numDataPoints) :
		AbstractDelayLineFunction(input, measuredOutput, numDataPoints) {
	}

	virtual void justifyParams(double *params) const {
		params[PARAM_L] = fmax(params[PARAM_L], minParamValue);
		params[PARAM_T] = fmax(params[PARAM_T], minParamValue);
	}

	// Creating a state-space representation using Rosenbrock system matrix
	virtual double getEstimatedValueAtPoint(int i, const double *params) const {
		// only positive values allowed (todo: choose the limits)
		double pL = fmax(params[PARAM_L], minParamValue);
		double pT = fmax(params[PARAM_T], minParamValue);

		// state-space params
		double lambda = exp(-1.0 / pT);

		// todo: find better initial value?
		double y = inputFunc->getValue(0) * params[PARAM_K];
		
		// The FO response function is indirect, so we need to iterate all previous values to find the current one
		for (int j = 0; j <= i; j++) {
			// delayed input
			double inp = inputFunc->getValue((double)j - pL);

			// indirect model response in Controllable Canonical Form (1st order CCF)
			y = lambda * y + params[PARAM_K] * (1.0 - lambda) * inp;
		}
		return y;
	}
};


// "Overdamped" SODPT indirect transfer function used for step response analytic simulation (xi > 1)
// The Laplace representation is: K * exp(-L * s) / ((T1*T2)*s^2 + (T1+T2)*s + 1)
class SecondOrderPlusDelayLineOverdampedFunction : public AbstractDelayLineFunction<4> {
public:
	SecondOrderPlusDelayLineOverdampedFunction(const InputFunction *input, const float *measuredOutput, int numDataPoints) :
		AbstractDelayLineFunction(input, measuredOutput, numDataPoints) {
	}

	virtual void justifyParams(double *params) const {
		params[PARAM_L] = fmax(params[PARAM_L], minParamValue);
		params[PARAM_T] = fmax(params[PARAM_T], minParamValue);
		params[PARAM_T2] = fmax(params[PARAM_T2], minParamValue);
	}

	// Creating a state-space representation using Rosenbrock system matrix
	virtual double getEstimatedValueAtPoint(int i, const double *params) const {
		// only positive values allowed (todo: choose the limits)
		double pL = fmax(params[PARAM_L], minParamValue);
		double pT = fmax(params[PARAM_T], minParamValue);
		double pT2 = fmax(params[PARAM_T2], minParamValue);

		// state-space params
		double lambda = exp(-1.0 / pT);
		double lambda2 = exp(-1.0 / pT2);

		// todo: find better initial values?
		double x = inputFunc->getValue(0) * params[PARAM_K];
		double y = inputFunc->getValue(0) * params[PARAM_K];

		// The SO response function is indirect, so we need to iterate all previous values to find the current one
		for (int j = 0; j <= i; j++) {
			// delayed input
			double inp = inputFunc->getValue((double)j - pL);

			// indirect model response in Controllable Canonical Form (2nd order CCF)
			y = lambda2 * y + (1.0 - lambda2) * x;
			x = lambda * x + params[PARAM_K] * (1.0 - lambda) * inp;
		}
		return y;
	}
};

