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

template <int numParams>
class LMSFunction {
public:
	// Get the total number of data points
	virtual double getNumPoints() const = 0;

	virtual void justifyParams(double *params) const = 0;

	/// Returns the y value of the function for the given x and vector of parameters
	virtual double getEstimatedValueAtPoint(int i, const double *params) const = 0;

	/// Returns the residual (error delta) of the function (return (dataPoints[i] - estimatedPoint(i)))
	virtual double getResidual(int i, const double *params) const = 0;

	/// Return the partial derivate of the function with respect to parameter pIndex at point [i].
	/// Can be overridden if analytical gradient function's representation is available
	virtual double getPartialDerivative(int i, const double *params, int pIndex) const {
		// some magic value
		const double delta = 1.0e-6;
		// we need to alter parameters around the neighborhood of 'pIndex', so we make a working copy
		double tmpParams[numParams];
		for (int k = 0; k < numParams; k++)
			tmpParams[k] = params[k];

		tmpParams[pIndex] = params[pIndex] + delta;
		double dplusResult = getEstimatedValueAtPoint(i, tmpParams);

		tmpParams[pIndex] = params[pIndex] - delta;
		double dminusResult = getEstimatedValueAtPoint(i, tmpParams);

		return (dplusResult - dminusResult) / (delta * 2.0);
	}
};

template<int numParams>
class LevenbergMarquardtSolver {
public:
	// ctor
	LevenbergMarquardtSolver(LMSFunction<numParams> *func, double parameters[numParams]) {
		this->func = func;
		this->parameters = parameters;
	}

	// lambda - magic coef.
	// maxIterations - if too many iterations (loop exit condition #1)
	// minDelta - if the progress on iteration is too small (loop exit condition #2)
	int solve(double lambda_ = 0.001, double minDelta = 1e-15, int maxIterations = 100) {
		this->lambda = lambda_;

		iterationCount = 0;
		
		double delta = 0;
		do {
			double merit = calcMerit(parameters);
			calcGradient();
			calcHessian();
			bool isSolved = calcNewParameters();
			double newMerit = calcMerit(newParameters);
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
			printf("[%d] (%g,%g,%g,%g) l=%g m=%g (%g-%g = %g)\r\n", iterationCount, parameters[0], parameters[1], parameters[2], parameters[3], lambda, merit,
				newMerit, merit, newMerit - merit);
#endif
			iterationCount++;
		} while (delta > minDelta && iterationCount < maxIterations);
		return iterationCount;
	}

	double *getParameters() const {
		return parameters;
	}
	
	// Calculate the sum of the squares of the residuals
	double calcMerit(double *params) {
		double res = 0;
		for (int i = 0; i < func->getNumPoints(); i++) {
			double r = func->getResidual(i, params);
			res += r * r;
		}
		return res;
	}

protected:
	// Find the parameter increments by solving the Hessian x Gradient equation
	bool calcNewParameters() {
		// get H^-1 matrix (inverse Hessian)
		double hinv[numParams][numParams];
		bool ret = MatrixHelper<double, numParams>::inverseMatrix(hinv, hessian);
		if (!ret)
			return false;

		for (int row = 0; row < numParams; row++) {
			double increment = 0;
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
				double res = 0;
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
			double res = 0;
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
	double *parameters; // [numParams]

	// Incremented (next step) parameters vector
	double newParameters[numParams];

	// Hessian matrix
	double hessian[numParams][numParams];
	// Gradients vector
	double gradient[numParams];

	// some magic number
	const double lambdaMultiplier = 10.0;

	// coeff used to adapt the descent step (and speed)
	double lambda;

	// Total number of iterations
	int iterationCount;
};
