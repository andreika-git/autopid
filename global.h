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


#ifndef CONTROLLERS_GENERATED_ENGINE_CONFIGURATION_GENERATED_STRUCTURES_H
struct pid_s {
	float pFactor;
	float iFactor;
	float dFactor;
	float offset;
	float periodMs;

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
