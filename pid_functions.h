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

// PID model params:
enum {
	PARAM_Kp = 0,
	PARAM_Ki,
	PARAM_Kd,
};
const int numParamsForPid = 3;

// the limit used for K and L params
static const double minParamKL = 0.00001;

// T (or T2) has a separate limit, because exp(-1/T) is zero for small T values 
// (and we need to compare between them to find a direction for optimization method).
static const double minParamT0 = 0.0001, minParamT1 = 0.01;

class InputFunction {
public:
	InputFunction(double timeScale_) : timeScale(timeScale_) {
	}

	// i = index, d = delay time (in seconds?)
	virtual double getValue(double i, double d) const = 0;

	virtual double getTimeScale() const {
		return timeScale;
	}

protected:
	// time scale needed to synchronize between the virtual step function and the real measured data
	//  timeScale=100 means 100 points per second.
	double timeScale;
};

// Heaviside step function interpolated between 'min' and 'max' values with 'stepPoint' time offset
class StepFunction : public InputFunction
{
public:
	StepFunction(double minValue_, double maxValue_, double offset_, double stepPoint_, double timeScale_) :
				minValue(minValue_), maxValue(maxValue_), offset(offset_), stepPoint(stepPoint_), InputFunction(timeScale_) {
	}

	virtual double getValue(double i, double d) const {
		// delayed index
		double id = i - d * timeScale;
#ifdef INTERPOLATED_STEP_FUNCTION
		// the delay parameter L may not be integer, so we have to interpolate between the closest input values (near and far in the past)
		int I = (int)id;
		double fract = id - I;	// 0 = choose near value, 1 = choose far value
		// find two closest input values for the given delay
		double vNear = (I < stepPoint) ? minValue : maxValue;
		double vFar = (I + 1 < stepPoint) ? minValue : maxValue;
		// interpolate
		return offset + vFar * fract + vNear * (1.0f - fract);
#else
		return offset + ((id < stepPoint) ? minValue : maxValue);
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
	StoredDataInputFunction(double timeScale_) : InputFunction(timeScale_) {
		reset();
	}

	void reset() {
		inputData.init();
	}

	void addDataPoint(float v) {
		// todo: support data scaling
		assert(inputData.getNumDataPoints() <= numPoints);
		inputData.addDataPoint(v);
	}

	virtual double getValue(double i, double d) const {
		return inputData.getValue((float)(i - d * timeScale));
	}

private:
	AveragingDataBuffer<numPoints> inputData;
};

// Abstract indirect transfer function used for step response analytic simulation
template <int numParams>
class AbstractDelayLineFunction : public LMSFunction<numParams> {
public:
	AbstractDelayLineFunction(const InputFunction *input, const float *measuredOutput, int numDataPoints, double minParamT) {
		dataPoints = measuredOutput;
		inputFunc = input;
		numPoints = numDataPoints;
		this->minParamT = minParamT;
	}

	virtual double getResidual(int i, const double *params) const {
		return dataPoints[i] - getEstimatedValueAtPoint(i, params);
	}

	virtual double getEstimatedValueAtPoint(int i, const double *params) const = 0;

	// Get the total number of data points
	virtual int getNumPoints() const {
		return numPoints;
	}

	float getDataPoint(int i) const {
		return dataPoints[i];
	}

protected:
	const InputFunction *inputFunc;
	const float *dataPoints;
	int numPoints;
	double minParamT;
};

// FODPT indirect transfer function used for step response analytic simulation.
// Used mostly as an approximate model for chemical or thermal processes
// The Laplace representation is: K * exp(-L*s) / (T*s + 1)
class FirstOrderPlusDelayLineFunction : public AbstractDelayLineFunction<3> {
public:
	FirstOrderPlusDelayLineFunction(const InputFunction *input, const float *measuredOutput, int numDataPoints, double minParamT) :
		AbstractDelayLineFunction(input, measuredOutput, numDataPoints, minParamT) {
	}

	virtual void justifyParams(double *params) const {
		params[PARAM_L] = fmax(params[PARAM_L], minParamKL);
		params[PARAM_T] = fmax(params[PARAM_T], minParamT);
		params[PARAM_K] = (fabs(params[PARAM_K]) < minParamKL) ? minParamKL : params[PARAM_K];
	}

	// Creating a state-space representation using Rosenbrock system matrix
	virtual double getEstimatedValueAtPoint(int i, const double *params) const {
		// only positive values allowed (todo: choose the limits)
		double pL = fmax(params[PARAM_L], minParamKL);
		double pT = fmax(params[PARAM_T], minParamT);
		double pK = (fabs(params[PARAM_K]) < minParamKL) ? minParamKL : params[PARAM_K];

		// state-space params
		double lambda = exp(-1.0 / (pT * inputFunc->getTimeScale()));

		// todo: find better initial value?
		double y = inputFunc->getValue(0, 0) * pK;
		
		// The FO response function is indirect, so we need to iterate all previous values to find the current one
		for (int j = 0; j <= i; j++) {
			// delayed input
			double inp = inputFunc->getValue((double)j, pL);

			// indirect model response in Controllable Canonical Form (1st order CCF)
			y = lambda * y + pK * (1.0 - lambda) * inp;
		}
		return y;
	}
};


// "Overdamped" SODPT indirect transfer function used for step response analytic simulation (xi > 1)
// Used mostly as an approximate model for electro-mechanical processes (e.g. manometer
// The Laplace representation is: K * exp(-L * s) / ((T1*T2)*s^2 + (T1+T2)*s + 1)
class SecondOrderPlusDelayLineOverdampedFunction : public AbstractDelayLineFunction<4> {
public:
	SecondOrderPlusDelayLineOverdampedFunction(const InputFunction *input, const float *measuredOutput, int numDataPoints, double minParamT) :
		AbstractDelayLineFunction(input, measuredOutput, numDataPoints, minParamT) {
	}

	virtual void justifyParams(double *params) const {
		params[PARAM_L] = fmax(params[PARAM_L], minParamKL);
		params[PARAM_T] = fmax(params[PARAM_T], minParamT);
		params[PARAM_T2] = fmax(params[PARAM_T2], minParamT);
		params[PARAM_K] = (fabs(params[PARAM_K]) < minParamKL) ? minParamKL : params[PARAM_K];
	}

	// Creating a state-space representation using Rosenbrock system matrix
	virtual double getEstimatedValueAtPoint(int i, const double *params) const {
		// only positive values allowed (todo: choose the limits)
		double pL = fmax(params[PARAM_L], minParamKL);
		double pT = fmax(params[PARAM_T], minParamT);
		double pT2 = fmax(params[PARAM_T2], minParamT);
		double pK = (fabs(params[PARAM_K]) < minParamKL) ? minParamKL : params[PARAM_K];

		// state-space params
		double lambda = exp(-1.0 / (pT * inputFunc->getTimeScale()));
		double lambda2 = exp(-1.0 / (pT2 * inputFunc->getTimeScale()));

		// todo: find better initial values?
		double x = inputFunc->getValue(0, 0) * pK;
		double y = inputFunc->getValue(0, 0) * pK;

		// The SO response function is indirect, so we need to iterate all previous values to find the current one
		for (int j = 0; j <= i; j++) {
			// delayed input
			double inp = inputFunc->getValue((double)j, pL);

			// indirect model response in Controllable Canonical Form (2nd order CCF)
			y = lambda2 * y + (1.0 - lambda2) * x;
			x = lambda * x + pK * (1.0 - lambda) * inp;
		}
		return y;
	}
};

// Harriot's relation function (based on the graph)
// Used to approximate initial parameters for SOPDT model
// See: findSecondOrderInitialParamsHarriot() and "Harriot P. Process control (1964). McGraw-Hill. USA."
class HarriotFunction {
public:
	double getValue(double x) const {
		return buf.getValue((float)((x - 2.8 / 719.0 - 0.26) * 719.0 / 2.8));
	}

	static constexpr double minX = 2.8 / 719.0 - 0.26;
	static constexpr double maxX = 33 / 719.0 * 2.8 + minX;

private:
	const AveragingDataBuffer<34> buf = { {
		0.500000000f, 0.589560440f, 0.624725275f, 0.652747253f, 0.675274725f, 0.694505495f, 0.712637363f, 0.729120879f,
		0.743406593f, 0.757142857f, 0.769780220f, 0.781318681f, 0.793956044f, 0.804395604f, 0.814285714f, 0.824725275f,
		0.834065934f, 0.844505495f, 0.853296703f, 0.862637363f, 0.870879121f, 0.880219780f, 0.889010989f, 0.897802198f,
		0.906593407f, 0.915384615f, 0.924175824f, 0.933516484f, 0.942857143f, 0.953296703f, 0.963736264f, 0.975274725f,
		0.986813187f, 1.000000000f
	}, 34, 0 };
};
