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
	GaussianFunction(int numPoints_, double_t *xValues_, double_t *yValues_) : numPoints(numPoints_), xValues(xValues_), yValues(yValues_) {
	}

	virtual void justifyParams(double_t *params) const {
	}

	// Get the total number of data points
	virtual int getNumPoints() const {
		return numPoints;
	}

	virtual double_t calcEstimatedValuesAtPoint(int pi, const double_t *params) const {
		double_t val = 0.0;
		for (int j = 0, i = 0; j < (numParams / 3); j++, i += 3)
		{
			double_t arg = (xValues[pi] - params[i + 1]) / params[i + 2];
			val += params[i] * exp(-arg * arg);
		}
		return val;
	}

	virtual double_t getResidual(int i, const double_t *params) const {
		return yValues[i] - getEstimatedValueAtPoint(i, params);
	}

private:
	int numPoints;
	double_t *xValues;
	double_t *yValues;
};

void testGaussianFunction() {
	int i;

	const int numParams = 6;
	const int numPoints = 100;
	
	const double_t goodParams[numParams] = { 5, 2, 3, 2, 5, 3 };
	
	double_t xValues[numPoints];
	for (i = 0; i < numPoints; i++) {
		xValues[i] = 0.1 * (double_t)(i + 1); 
	}
	
	double_t yValues[numPoints];

	GaussianFunction func(numPoints, xValues, yValues);

	// imitate "real" data by using our ideal function
	for (i = 0; i < numPoints; i++) {
		yValues[i] = func.getEstimatedValueAtPoint(i, goodParams);
	}

	double_t params[numParams] = { 4, 2, 2, 2, 5, 2 };
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
	double_t a[n][n];
	// fill with some arbitrary values
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j] = 1.0 / (i * n + j + 1);
		}
	}

	double_t ai[n][n];
	bool ret = MatrixHelper<double_t, n>::inverseMatrix(ai, a);
	ASSERT_EQ(true, ret);

	double_t mul[n][n];
	MatrixHelper<double_t, n>::multiplyMatrix(mul, ai, a);
	// A^-1 * A = I
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			double_t va = (i == j) ? 1.0 : 0;	// identity matrix element
			ASSERT_DOUBLE_EQ(va, mul[i][j]);
		}
	}
}

