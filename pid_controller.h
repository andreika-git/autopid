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

	double_t limitOutput(double_t v) {
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

	double_t getOutput(double_t target, double_t input, double_t dTime) {
		double_t error = target - input;
		pTerm = p.pFactor * error;
		iTerm += p.iFactor * dTime * error;
		dTerm = p.dFactor / dTime * (error - previousError);
		previousError = error;

		return limitOutput(pTerm + iTerm + dTerm + p.offset);
	}

protected:
	double_t pTerm, iTerm, dTerm;
	double_t previousError = 0;
};

// PID with derivative filtering (backward differences) and integrator anti-windup.
// See: Wittenmark B., Åström K., Årzén K. IFAC Professional Brief. Computer Control: An Overview. 
// Two additional parameters used: derivativeFilterLoss and antiwindupFreq
//   (If both are 0, then this controller is identical to PidParallelController)
class PidIndustrialController : public PidParallelController {
public:
	PidIndustrialController(const pid_s & p_) : PidParallelController(p_) {
	}

	double_t getOutput(double_t target, double_t input, double_t dTime) {
		double_t ad, bd;
		double_t error = target - input;
		pTerm = p.pFactor * error;

		// update the I-term
		iTerm += p.iFactor * dTime * error;
		
		// calculate dTerm coefficients
		if (fabs(p.derivativeFilterLoss) > DBL_EPSILON) {
			// restore Td in the Standard form from the Parallel form: Td = Kd / Kc
			double_t Td = p.dFactor / p.pFactor;
			// calculate the backward differences approximation of the derivative term
			ad = Td / (Td + dTime / p.derivativeFilterLoss);
			bd = p.pFactor * ad / p.derivativeFilterLoss;
		} else {
			// According to the Theory of limits, if p.derivativeFilterLoss -> 0, then 
			//   lim(ad) = 0; lim(bd) = p.pFactor * Td / dTime = p.dFactor / dTime
			//   i.e. dTerm becomes equal to PidParallelController's
			ad = 0.0;
			bd = p.dFactor / dTime;
		}
		
		// (error - previousError) = (target-input) - (target-prevousInput) = -(input - prevousInput)
		dTerm = dTerm * ad + (error - previousError) * bd;

		// calculate output and apply the limits
		double_t output = pTerm + iTerm + dTerm + p.offset;
		double_t limitedOutput = limitOutput(output);

		// apply the integrator anti-windup
		// If p.antiwindupFreq = 0, then iTerm is equal to PidParallelController's
		iTerm += dTime * p.antiwindupFreq * (limitedOutput - output);
		
		// update the state
		previousError = error;

		return limitedOutput;
	}
};

// C(s) = Kp + (Ki / s) + (N * Kd * s / (1 + N / s))
// The Integral term is discretized using backward Euler method
// See: https://www.scilab.org/discrete-time-pid-controller-implementation
class PidDerivativeFilterController : public PidController {
public:
	PidDerivativeFilterController(const pid_s & p_, double_t n_) : PidController(p_), N(n_) {
	}

	double_t getOutput(double_t target, double_t input, double_t dTime) {
		double_t error = target - input;
		double_t a0 = (1.0 + N * dTime);
		double_t a1 = -(2.0 + N * dTime);

		double_t a2 = 1.0;
		double_t b0 = p.pFactor * (1.0 + N * dTime) + p.iFactor * dTime * (1.0 + N * dTime) + p.dFactor * N;
		double_t b1 = -(p.pFactor * (2.0 + N * dTime) + p.iFactor * dTime + 2.0 * p.dFactor * N);
		double_t b2 = p.pFactor + p.dFactor * N;

		double_t ku1 = a1 / a0; 
		double_t ku2 = a2 / a0; 
		double_t ke0 = b0 / a0; 
		double_t ke1 = b1 / a0; 
		double_t ke2 = b2 / a0;

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
	double_t e2 = 0, e1 = 0, e0 = 0, u2 = 0, u1 = 0, u0 = 0;
	double_t N = 1;
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

	void addPoint(double_t i, double_t value, double_t target) {
		double_t e = target - value;
		itae += i * fabs(e);
		ise += e * e;
		double_t overshoot = (value - target) / target;
		if (overshoot > 0 && overshoot > maxOvershoot)
			maxOvershoot = overshoot;
		lastValue = value;
	}

	double_t getItae() const {
		return itae;
	}

	double_t getIse() const {
		return ise;
	}

	double_t getMaxOvershoot() const {
		return maxOvershoot;
	}

	double_t getLastValue() const {
		return lastValue;
	}

	double_t getMerit() const {
		//return getItae();
#if 1
		double overShootWeight = 10000.0;
		return getItae() + overShootWeight * getMaxOvershoot() * getMaxOvershoot();
#endif
	}

private:
	double_t itae = 0;	// Integral time-weighted absolute error
	double_t ise = 0;		// Integral square error
	double_t maxOvershoot = 0;
	double_t lastValue = 0;
};