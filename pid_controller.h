/*
* @file	pid_controller.h
*
* PID Controller models needed to verify the parameters.
*
* @date Oct 02, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"

class PidController {
public:
	PidController(const pid_s & p_) : p(p_) {
	}

	double limitOutput(double v) {
		if (v < p.minValue)
			v = p.minValue;
		if (v > p.maxValue)
			v = p.maxValue;
		return v;
	}

protected:
	const pid_s p;
};

class PidParallelController : public PidController {
public:
	PidParallelController(const pid_s & p_) : PidController(p_) {
		pTerm = iTerm = dTerm = 0.0;
	}

	double getOutput(double target, double input, double dTime) {
		double error = target - input;
		pTerm = p.pFactor * error;
		iTerm += p.iFactor * dTime * error;
		dTerm = p.dFactor / dTime * (error - previousError);
		previousError = error;

		return limitOutput(pTerm + iTerm + dTerm + p.offset);
	}

protected:
	double pTerm, iTerm, dTerm;
	double previousError = 0;
};

// C(s) = Kp + (Ki / s) + (N * Kd * s / (1 + N / s))
// The Integral term is discretized using backward Euler method
// See: https://www.scilab.org/discrete-time-pid-controller-implementation
class PidDerivativeFilterController : public PidController {
public:
	PidDerivativeFilterController(const pid_s & p_, double n_) : PidController(p_), N(n_) {
	}

	double getOutput(double target, double input, double dTime) {
		double error = target - input;
		double a0 = (1.0 + N * dTime);
		double a1 = -(2.0 + N * dTime);

		double a2 = 1.0;
		double b0 = p.pFactor * (1.0 + N * dTime) + p.iFactor * dTime * (1.0 + N * dTime) + p.dFactor * N;
		double b1 = -(p.pFactor * (2.0 + N * dTime) + p.iFactor * dTime + 2.0 * p.dFactor * N);
		double b2 = p.pFactor + p.dFactor * N;

		double ku1 = a1 / a0; 
		double ku2 = a2 / a0; 
		double ke0 = b0 / a0; 
		double ke1 = b1 / a0; 
		double ke2 = b2 / a0;

		e2 = e1; 
		e1 = e0; 
		u2 = u1;
		u1 = u0;

		e0 = error;
		u0 = -ku1 * u1 - ku2 * u2 + ke0 * e0 + ke1 * e1 + ke2 * e2;

		u0 = limitOutput(u0);
		
		return u0;
	}

protected:
	double e2 = 0, e1 = 0, e0 = 0, u2 = 0, u1 = 0, u0 = 0;
	double N = 1;
};

// Calculate ITAE/ISE and Overshoot
class PidAccuracyMetric {
public:
	void reset() {
		itae = 0;
		ise = 0;
		maxOvershoot = 0;
		lastValue = 0;
	}

	void addPoint(double i, double value, double target) {
		double e = target - value;
		itae += i * fabs(e);
		ise += e * e;
		double overshoot = (value - target) / target;
		if (overshoot > 0 && overshoot > maxOvershoot)
			maxOvershoot = overshoot;
		lastValue = value;
	}

	double getItae() const {
		return itae;
	}

	double getIse() const {
		return ise;
	}

	double getMaxOvershoot() const {
		return maxOvershoot;
	}

	double getLastValue() const {
		return lastValue;
	}

private:
	double itae = 0;	// Integral time-weighted absolute error
	double ise = 0;		// Integral square error
	double maxOvershoot = 0;
	double lastValue = 0;
};