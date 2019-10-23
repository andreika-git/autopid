// Dummy file replacement

#pragma once


#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gtest/gtest.h"
#include "gmock/gmock.h"


//!!!!!!!!!!!!!!
#define float_t double_t

#ifndef CONTROLLERS_GENERATED_ENGINE_CONFIGURATION_GENERATED_STRUCTURES_H
struct pid_s {
	float_t pFactor;
	float_t iFactor;
	float_t dFactor;
	float_t offset;
	float_t periodMs;

	float_t antiwindupFreq;		// = 1/ResetTime
	float_t derivativeFilterLoss;	// = 1/Gain

	int16_t minValue;
	int16_t maxValue;
};
typedef struct pid_s pid_s;
#endif

#if 1	// remove!
#define LMS_DEBUG
#define PID_DEBUG
#endif

#define EPS1D 0.1
#define EPS2D 0.01
#define EPS3D 0.001
#define EPS4D 0.0001
#define EPS5D 0.00001
