/*
* @file	matrix_helper.h
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"


template<typename FLOAT, int N>
class MatrixHelper {
public:
	// Gauss-Jordan elimination method
	static bool inverseMatrix(FLOAT dst[N][N], FLOAT src[N][N]) {
		// a working copy of src (needs modification)
		FLOAT tmp[N][N];

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				tmp[i][j] = src[i][j];
				dst[i][j] = (i == j) ? 1.0 : 0;	// fill as an identity matrix first
			}
		}

		// determinant
		FLOAT det = 1.0;

		// for each pass, find the maximum element in the pivot column.
		for (int pass = 0; pass < N; pass++) {
			int row, maxRow = pass;
			// find the largest element in the column
			for (row = pass; row < N; row++) {
				if (fabs(tmp[row][pass]) > fabs(tmp[maxRow][pass]))
					maxRow = row;
			}
			// interchange the elements of these rows
			if (maxRow != pass) {
				for (int col = 0; col < N; col++) {
					FLOAT temp = dst[pass][col];
					dst[pass][col] = dst[maxRow][col];
					dst[maxRow][col] = temp;

					if (col >= pass) {
						temp = tmp[pass][col];
						tmp[pass][col] = tmp[maxRow][col];
						tmp[maxRow][col] = temp;
					}
				}
			}

			// calculate the determinant as the product of the elements
			FLOAT coef = tmp[pass][pass];
			if (coef == 0.0)
				return false;
			det *= coef;

			// normalize
			for (int col = 0; col < N; col++) {
				dst[pass][col] /= coef;
				if (col >= pass)
					tmp[pass][col] /= coef;
			}

			// add a multiple of the pivot row to each row
			for (row = 0; row < N; row++) {
				if (row != pass) {
					coef = tmp[row][pass];
					for (int col = 0; col < N; col++) {
						dst[row][col] -= coef * dst[pass][col];
						tmp[row][col] -= coef * tmp[pass][col];
					}
				}
			}
		}

		return true;
	}

	static void multiplyMatrix(FLOAT dst[N][N], FLOAT src1[N][N], FLOAT src2[N][N]) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				dst[i][j] = 0;
				for (int k = 0; k < N; k++) {
					dst[i][j] += src1[i][k] * src2[k][j];
				}
			}
		}
	}

};


