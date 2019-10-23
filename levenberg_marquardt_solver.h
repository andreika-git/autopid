/*
* @file	levenberg_marquardt_solver.h
*
* Levenberg-Marquardt Algorithm, and efficient non-linear optimization solver used for regression analysis ("least squares problem").
* It basically combines Gauss-Newton method and gradient descent, but using an approximation for computing a Hessian matrix!
* The code is based on "W.Press, S.Teukolsky, W.Vetterling, B.Flannery. Numerical Recipes in Fortran 77, 1992."
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"

#include "matrix_helper.h"

#ifndef LM_USE_CACHE
#define LM_USE_CACHE
#endif /* LM_USE_CACHE */

// Used to get all the values for the Levenberg-Marquardt solver.
template <int numParams>
class LMSFunction {
public:
	// Get the total number of data points
	virtual int getNumPoints() const = 0;

	virtual void justifyParams(double_t *params) const = 0;

	/// Returns the y value of the function for the given x and vector of parameters
	virtual double_t getEstimatedValueAtPoint(int i, const double_t *params) const {
#ifdef LM_USE_CACHE
		return points[0][i];
#else
		return calcEstimatedValuesAtPoint(i, params);
#endif
	}

	/// Returns the y value of the function for the given x and vector of parameters
	virtual double_t calcEstimatedValuesAtPoint(int i, const double_t *params) const = 0;

	/// Returns the residual (error delta) of the function (return (dataPoints[i] - estimatedPoint(i)))
	virtual double_t getResidual(int i, const double_t *params) const = 0;

	// Calculate the sum of the squares of the residuals
	virtual double_t calcMerit(double_t *params) const {
		double_t res = 0;
		for (int i = 0; i < getNumPoints(); i++) {
			double_t r = getResidual(i, params);
			res += r * r;
		}
		return res;
	}

	/// Return the partial derivate of the function with respect to parameter pIndex at point [i].
	/// Can be overridden if analytical gradient function's representation is available
	virtual double_t getPartialDerivative(int i, const double_t *params, int pIndex) const {
		// we need to alter parameters around the neighborhood of 'pIndex', so we make a working copy
#ifdef LM_USE_CACHE
		int p = 1 + pIndex * 2;
		double_t dplusResult = points[p][i];
		double_t dminusResult = points[p + 1][i];
#else
		double_t tmpParams[numParams];
		for (int k = 0; k < numParams; k++)
			tmpParams[k] = params[k];

		tmpParams[pIndex] = params[pIndex] + delta;
		double_t dplusResult = getEstimatedValueAtPoint(i, tmpParams);

		tmpParams[pIndex] = params[pIndex] - delta;
		double_t dminusResult = getEstimatedValueAtPoint(i, tmpParams);
#endif

		return (dplusResult - dminusResult) / (delta * 2.0);
	}

	void calculateAllPoints(const double_t *params) {
#ifdef LM_USE_CACHE
		for (int i = 0; i < numPointGroups; i++) {
			points[i].resize(getNumPoints());
		}
		// calculate displaced points for partial derivatives
		for (int pIndex = 0; pIndex < numParams; pIndex++) {
			double_t tmpParams[numParams];
			for (int k = 0; k < numParams; k++)
				tmpParams[k] = params[k];
			// 2 points for each derivative: +delta and -delta
			for (int plusMinus = 0; plusMinus < 2; plusMinus++) {
				tmpParams[pIndex] = params[pIndex] + (plusMinus == 0 ? delta : -delta);
				int p = 1 + pIndex * 2 + plusMinus;
				for (int i = 0; i < getNumPoints(); i++) {
					points[p][i] = calcEstimatedValuesAtPoint(i, tmpParams);
				}
			}
			for (int i = 0; i < getNumPoints(); i++) {
				points[0][i] = calcEstimatedValuesAtPoint(i, params);
			}
		}
#endif
	}

private:
#ifdef LM_USE_CACHE
	static const int numPointGroups = numParams * 2 + 1;
	std::vector<double_t> points[numPointGroups];
#endif

protected:
	// some magic value for the differentiation step to calculate a partial derivative
	double_t delta = 1.0e-6;
};

template<int numParams>
class LevenbergMarquardtSolver {
public:
	// ctor
	LevenbergMarquardtSolver(LMSFunction<numParams> *func, double_t *parameters) {
		this->func = func;
		this->parameters = parameters;
	}

	// lambda - magic coef.
	// maxIterations - if too many iterations (loop exit condition #1)
	// minDelta - if the progress on iteration is too small (loop exit condition #2)
	int solve(double_t lambda_ = 0.001, double_t minDelta = 1e-15, int maxIterations = 100) {
		this->lambda = lambda_;

		iterationCount = 0;
		
		func->calculateAllPoints(parameters);
		double_t delta = 0;
		do {
			double_t merit = func->calcMerit(parameters);

			calcGradient();
			calcHessian();
			bool isSolved = calcNewParameters();
			func->calculateAllPoints(newParameters);
			double_t newMerit = func->calcMerit(newParameters);
			if (!isSolved) {
				return -iterationCount;
			}
			// if we don't like the new parameters
			if (newMerit >= merit) {
				// decrease the step
				lambda *= lambdaMultiplier;
			}
			// if we're good, accept them
			else {
				// update via copy
				memcpy(parameters, newParameters, sizeof(newParameters));
				// let's increase the step even more
				lambda /= lambdaMultiplier;
			}
			// find out if we progressed enough in this iteration
			delta = fabs(newMerit - merit);
#ifdef LMS_DEBUG
			printf("[%d] (%Lg,%Lg,%Lg,%Lg) l=%Lg merit %Lg->%Lg, dm=%Lg\r\n", iterationCount, (long double)parameters[0], (long double)parameters[1], (long double)parameters[2], (long double)parameters[3], 
				(long double)lambda, (long double)merit, (long double)newMerit, (long double)(newMerit - merit));
#endif
			iterationCount++;
		} while (delta > minDelta && iterationCount < maxIterations);
		return iterationCount;
	}

	double_t *getParameters() const {
		return parameters;
	}
	
protected:
	// Find the parameter increments by solving the Hessian x Gradient equation
	bool calcNewParameters() {
		// get H^-1 matrix (inverse Hessian)
		double_t hinv[numParams][numParams];
		bool ret = MatrixHelper<double_t, numParams>::inverseMatrix(hinv, hessian);
		if (!ret)
			return false;

		for (int row = 0; row < numParams; row++) {
			double_t increment = 0;
			for (int col = 0; col < numParams; col++) {
				increment += hinv[row][col] * gradient[col];
			}

			newParameters[row] = parameters[row] + increment;
		}
		func->justifyParams(newParameters);
		return true;
	}

	/// Calculate the Hessian matrix (2nd derivative) approximation
	void calcHessian() {
		for (int row = 0; row < numParams; row++) {
			for (int col = 0; col < numParams; col++) {
				double_t res = 0;
				for (int i = 0; i < func->getNumPoints(); i++) {
					res += func->getPartialDerivative(i, parameters, row) * func->getPartialDerivative(i, parameters, col);
				}
				hessian[row][col] = (row == col) ? res * (lambda + 1.0) : res;
			}
		}
	}

	// Calculate the 1st derivatives of the residual func
	void calcGradient() {
		for (int row = 0; row < numParams; row++) {
			double_t res = 0;
			for (int i = 0; i < func->getNumPoints(); i++) {
				res += func->getResidual(i, parameters) * func->getPartialDerivative(i, parameters, row);
			}
			gradient[row] = res;
		}
	}

private:
	// optimization function
	LMSFunction<numParams> *func;

	// Current (accepted) parameters vector
	double_t *parameters; // [numParams]

	// Incremented (next step) parameters vector
	double_t newParameters[numParams];

	// Hessian matrix
	double_t hessian[numParams][numParams];
	// Gradients vector
	double_t gradient[numParams];

	// some magic number
	const double_t lambdaMultiplier = 10.0;

	// coeff used to adapt the descent step (and speed)
	double_t lambda;

	// Total number of iterations
	int iterationCount;
};
