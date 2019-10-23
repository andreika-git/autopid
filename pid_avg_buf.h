/*
* @file	pid_avg_buf.h
*
* Chien-Hrones-Reswick Algorithm, a PID coefficient finder method using step-response measured data
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#pragma once

#include "global.h"

// Used to store measured data in the memory-limited buffer.
// The buffer adapts to the data size automatically by averaging stored values.
template<int maxPoints>
class AveragingDataBuffer {
public:
	// no default ctors!

	void init() {
		// zero buffer
		memset(buf, 0, sizeof(buf));
		num = 0;
		scaleShift = 0;
	}

	void addDataPoint(float_t v) {
		int idx;

		for (;;) {
			idx = num >> scaleShift;

			if (idx < maxPoints)
				break;
			// we're here because the buffer size is too small to hold the new point
			int idxHalf = idx / 2;
			// shrink the buffer twice using averaging, and clear the rest of the buffer
			for (int i = 0; i < maxPoints; i++) {
				buf[i] = (i < idxHalf) ? (buf[i * 2] + buf[i * 2 + 1]) * 0.5f : 0;
			}
			scaleShift++;
		}
		float_t numInTheCell = (float_t)(num - ((num >> scaleShift) << scaleShift));
		buf[idx] *= numInTheCell;
		buf[idx] += v;
		buf[idx] /= (numInTheCell + 1.0f);
		num++;
	}

	int getNumDataPoints() const {
		return (num < 1) ? 1 : ((num - 1) >> scaleShift) + 1;
	}

	float_t const *getBuf() const {
		return buf;
	}

	// This is a robust method for all rational indices
	float_t getValue(float_t i) const {
		// todo: this works only for scale=1 :(
		assert(scaleShift == 0);

		// reject empty buffer
		if (num < 1)
			return 0.0f;
		// singular?
		if (num == 1)
			return buf[0];

		int I = (int)i;
		// extrapolate to the left?
		if (I < 0)
			return buf[0];
		// extrapolate to the right?
		if (I >= num - 1)
			return buf[num - 1];
		// get 2 closest values and interpolate
		float_t fract = i - I;
		return buf[I + 1] * fract + buf[I] * (1.0f - fract);
	}

	float_t getAveragedData(int from, int to) const {
		float_t avg = 0.0f;
		for (int i = from; i <= to; i++) {
			avg += buf[i];
		}
		avg /= (float_t)(to - from + 1);
		return avg;
	}

	int findDataAt(float_t v, int from, int to) const {
		for (int i = from; i <= to; i++) {
			if (buf[i] > v)
				return i;
		}
		return -1;
	}

	// integrating using simple trapezoidal method
	float_t getArea(int from, int to, float_t v0) const {
		float_t sum = 0.0f;
		for (int i = from; i < to; i++) {
			sum += buf[i] + buf[i + 1] - 2.0f * v0;
		}
		// assume the time step is 1.0
		return sum * 0.5f;
	}

public:
	float_t buf[maxPoints];
	int num = 0;
	int scaleShift = 0;
};

