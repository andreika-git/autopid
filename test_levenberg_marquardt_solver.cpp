/*
* @file	test_levenberg_marquardt_solver.cpp
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#include "global.h"

#include "matrix_helper.h"
#include "levenberg_marquardt_solver.h"

// Use simple and well-known Gaussian function to test the solver
class GaussianFunction : public LMSFunction<6> {
	const int numParams = 6;
public:
	GaussianFunction(int numPoints_, double *xValues_, double *yValues_) : numPoints(numPoints_), xValues(xValues_), yValues(yValues_) {
	}

	virtual void justifyParams(double *params) const {
	}

	// Get the total number of data points
	virtual int getNumPoints() const {
		return numPoints;
	}

	virtual double getEstimatedValueAtPoint(int pi, const double *params) const {
		double val = 0.0;
		for (int j = 0, i = 0; j < (numParams / 3); j++, i += 3)
		{
			double arg = (xValues[pi] - params[i + 1]) / params[i + 2];
			val += params[i] * exp(-arg * arg);
		}
		return val;
	}

	virtual double getResidual(int i, const double *params) const {
		return yValues[i] - getEstimatedValueAtPoint(i, params);
	}

private:
	int numPoints;
	double *xValues;
	double *yValues;
};

void testGaussianFunction() {
	int i;

	const int numParams = 6;
	const int numPoints = 100;
	
	const double goodParams[numParams] = { 5, 2, 3, 2, 5, 3 };
	
	double xValues[numPoints];
	for (i = 0; i < numPoints; i++) {
		xValues[i] = 0.1 * (double)(i + 1); 
	}
	
	double yValues[numPoints];

	GaussianFunction func(numPoints, xValues, yValues);

	// imitate "real" data by using our ideal function
	for (i = 0; i < numPoints; i++) {
		yValues[i] = func.getEstimatedValueAtPoint(i, goodParams);
	}

	double params[numParams] = { 4, 2, 2, 2, 5, 2 };
	LevenbergMarquardtSolver<numParams> solver(&func, params);

	int iterationCount = solver.solve(0.001, 1e-15, 100);

	for (i = 0; i < numParams; i++) {
		ASSERT_DOUBLE_EQ(goodParams[i], params[i]);
	}

	ASSERT_EQ(9, iterationCount);
}

void testMatrixInverse() {
	int i, j;
	const int n = 4;
	double a[n][n];
	// fill with some arbitrary values
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j] = 1.0 / (i * n + j + 1);
		}
	}

	double ai[n][n];
	bool ret = MatrixHelper<double, n>::inverseMatrix(ai, a);
	ASSERT_EQ(true, ret);

	double mul[n][n];
	MatrixHelper<double, n>::multiplyMatrix(mul, ai, a);
	// A^-1 * A = I
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			double va = (i == j) ? 1.0 : 0;	// identity matrix element
			ASSERT_DOUBLE_EQ(va, mul[i][j]);
		}
	}
}

