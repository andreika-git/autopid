/*
* @file	test_chs_pid.cpp
*
* @date Sep 27, 2019
* @author andreika, (c) 2019
*/

#include "global.h"

#include "pid_open_loop_models.h"
#include "pid_controller.h"

extern void testGaussianFunction();
extern void testMatrixInverse();


TEST(pidAutoTune, testMeasuredDataBuffer) {
	const int numPoints = 2;
	AveragingDataBuffer<numPoints> buf;
	buf.init();
	for (int i = 0; i < 16; i++) {
		buf.addDataPoint((float)(i + 1));
		float v = (float)(1 << (int)(log(i) / log(2.0)));
		float v1 = 0.5f * (1.0f + v);
		float v2 = (i < 1) ? 0 : v1 * 2.0f + (i - v) * 0.5f;
		ASSERT_EQ(v1, buf.getBuf()[0]);
		ASSERT_EQ(v2, buf.getBuf()[1]);
		ASSERT_EQ((i == 0) ? 1 : 2, buf.getNumDataPoints());
	}
}

TEST(pidAutoTune, testFOPDT) {
	StepFunction stepFunc(0, 100, 0, 10, 1.0);
	FirstOrderPlusDelayLineFunction func(&stepFunc, nullptr, 0);
	double params[3];
	params[PARAM_K] = 2.0;
	params[PARAM_T] = 3.0;
	params[PARAM_L] = 4.0;

	double v = func.getEstimatedValueAtPoint(24, params);
	ASSERT_DOUBLE_EQ(25.200031003972988, v);
}

TEST(pidAutoTune, testSOPDTOverdamped) {
	StepFunction stepFunc(0, 100, 0, 10, 1.0);
	SecondOrderPlusDelayLineOverdampedFunction func(&stepFunc, nullptr, 0);
	double params[4];
	params[PARAM_K] = 2.0;
	params[PARAM_T] = 3.0;
	params[PARAM_T2] = 0.3;
	params[PARAM_L] = 4.0;

	double v = func.getEstimatedValueAtPoint(24, params);
	ASSERT_DOUBLE_EQ(25.200031003972988, v);
}


static const float outputData[] = { 13.29, 13.29, 13.33, 13.33, 13.33, 13.33, 13.33, 13.22, 13.22, 13.22, 13.22, 13.3, 13.3, 13.3, 13.3, 13.3, 13.34, 13.34, 13.34, 13.34, 13.34, 13.2, 13.2, 13.2, 13.2, 13.2, 13.29, 13.29, 13.29, 13.29, 13.29, 13.32, 13.32, 13.32, 13.32, 13.32, 13.19, 13.19, 13.19, 13.19, 13.19, 13.28, 13.28, 13.28, 13.28, 13.28, 13.32, 13.32, 13.32, 13.32, 13.32, 13.18, 13.18, 13.18, 13.18, 13.18, 13.27, 13.27, 13.27, 13.27, 13.27, 13.32, 13.32, 13.32, 13.32, 13.17, 13.17, 13.17, 13.17, 13.17, 13.27, 13.27, 13.27, 13.27, 13.27, 13.32, 13.32, 13.32, 13.32, 13.32, 13.16, 13.16, 13.16, 13.16, 13.16, 13.25, 13.25, 13.25, 13.25, 13.25, 13.3, 13.3, 13.3, 13.3, 13.3, 13.14, 13.14, 13.14, 13.14, 13.23, 13.23, 13.23, 13.23, 13.23, 13.28, 13.28, 13.28, 13.28, 13.28, 13.16, 13.16, 13.16, 13.16, 13.16, 13.25, 13.25, 13.25, 13.25, 13.25, 13.28, 13.28, 13.28, 13.28, 13.28, 13.2, 13.2, 13.2, 13.2, 13.2, 13.27, 13.27, 13.27, 13.27, 13.27, 13.29, 13.29, 13.29, 13.29, 13.24, 13.24, 13.24, 13.24, 13.24, 13.3, 13.3, 13.3, 13.3, 13.3, 13.3, 13.3, 13.3, 13.3, 13.29, 13.29, 13.29, 13.29, 13.29, 13.33, 13.33, 13.33, 13.33, 13.33, 13.3, 13.3, 13.3, 13.3, 13.3, 13.33, 13.33, 13.33, 13.33, 13.33, 13.36, 13.36, 13.36, 13.36, 13.36, 13.31, 13.37, 13.37, 13.37, 13.44, 13.44, 13.44, 13.44, 13.44, 13.44, 13.44, 13.44, 13.44, 13.44, 13.54, 13.54, 13.54, 13.54, 13.54, 13.62, 13.62, 13.62, 13.62, 13.62, 13.56, 13.56, 13.56, 13.56, 13.56, 13.68, 13.68, 13.68, 13.68, 13.76, 13.76, 13.76, 13.76, 13.76, 13.65, 13.65, 13.65, 13.65, 13.65, 13.78, 13.78, 13.78, 13.78, 13.78, 13.84, 13.84, 13.84, 13.84, 13.84, 13.95, 13.95, 13.95, 13.95, 13.95, 14.04, 14.04, 14.04, 14.04, 14.04, 13.91, 13.91, 13.91, 13.91, 13.91, 14.06, 14.06, 14.06, 14.06, 14.06, 14.11, 14.11, 14.11, 14.11, 14.23, 14.23, 14.23, 14.23, 14.23, 14.33, 14.33, 14.33, 14.33, 14.33, 14.37, 14.37, 14.37, 14.37, 14.37, 14.48, 14.48, 14.48, 14.48, 14.48, 14.36, 14.36, 14.36, 14.36, 14.36, 14.53, 14.53, 14.53, 14.53, 14.53, 14.59, 14.59, 14.59, 14.59, 14.59, 14.74, 14.74, 14.74, 14.74, 14.74, 14.85, 14.85, 14.85, 14.85, 14.85, 14.94, 14.94, 14.94, 14.94, 15.05, 15.05, 15.05, 15.05, 15.05, 14.91, 14.91, 14.91, 14.91, 14.91, 15.06, 15.06, 15.06, 15.06, 15.06, 15.05, 15.05, 15.05, 15.05, 15.05, 15.18, 15.18, 15.18, 15.18, 15.18, 15.23, 15.23, 15.23, 15.23, 15.23, 15.34, 15.34, 15.34, 15.34, 15.34, 15.4, 15.4, 15.4, 15.4, 15.4, 15.42, 15.42, 15.42, 15.42, 15.49, 15.49, 15.49, 15.49, 15.49, 15.32, 15.32, 15.32, 15.32, 15.32, 15.45, 15.45, 15.45, 15.45, 15.45, 15.43, 15.43, 15.43, 15.43, 15.43, 15.53, 15.53, 15.53, 15.53, 15.53, 15.58, 15.58, 15.58, 15.58, 15.58, 15.63, 15.63, 15.63, 15.63, 15.63, 15.67, 15.67, 15.67, 15.67, 15.67, 15.5, 15.5, 15.5, 15.5, 15.5, 15.61, 15.61, 15.61, 15.61, 15.61, 15.57, 15.57, 15.57, 15.57, 15.57, 15.66, 15.66, 15.66, 15.66, 15.66, 15.7, 15.7, 15.7, 15.7, 15.74, 15.74, 15.74, 15.74, 15.74, 15.77, 15.77, 15.77, 15.77, 15.77, 15.63, 15.63, 15.63, 15.63, 15.63, 15.7, 15.7, 15.7, 15.7, 15.59, 15.59, 15.59, 15.59, 15.68, 15.68, 15.68, 15.68, 15.68, 15.68, 15.68, 15.68, 15.68, 15.68, 15.75, 15.75, 15.75, 15.75, 15.75, 15.77, 15.77, 15.77, 15.77, 15.77, 15.8, 15.8, 15.8, 15.8, 15.8, 15.83, 15.83, 15.83, 15.83, 15.83, 15.71, 15.71, 15.71, 15.71, 15.77, 15.77, 15.77, 15.77, 15.77, 15.57, 15.57, 15.57, 15.57, 15.57, 15.68, 15.68, 15.68, 15.68, 15.68, 15.61, 15.61, 15.61, 15.61, 15.61, 15.71, 15.71, 15.71, 15.71, 15.71, 15.7, 15.7, 15.7, 15.7, 15.7, 15.77, 15.77, 15.77, 15.77, 15.79, 15.79, 15.79, 15.79, 15.79, 15.82, 15.82, 15.82, 15.82, 15.82, 15.84, 15.84, 15.84, 15.78, 15.78, 15.78, 15.78, 15.83, 15.83, 15.83, 15.83, 15.83, 15.62, 15.62, 15.62, 15.62, 15.73, 15.73, 15.73, 15.73, 15.73, 15.66, 15.66, 15.66, 15.66, 15.66, 15.75, 15.75, 15.75, 15.75, 15.75, 15.75, 15.75, 15.75, 15.75, 15.75, 15.81, 15.81, 15.81, 15.81, 15.81, 15.84, 15.84, 15.84, 15.84, 15.84, 15.86, 15.86, 15.86, 15.86, 15.86, 15.88, 15.88, 15.88, 15.88, 15.88, 15.74, 15.74, 15.74, 15.74, 15.74, 15.81, 15.81, 15.81, 15.81, 15.81, 15.69, 15.69, 15.69, 15.69, 15.69, 15.77, 15.77, 15.77, 15.77, 15.77, 15.79, 15.79, 15.79, 15.79, 15.86, 15.86, 15.86, 15.86, 15.86, 15.88, 15.88, 15.88, 15.88, 15.88, 15.82, 15.82, 15.82, 15.82, 15.82, 15.88, 15.88, 15.88, 15.88, 15.88, 15.69, 15.69, 15.69, 15.69, 15.69, 15.78, 15.78, 15.78, 15.78, 15.78, 15.79, 15.79, 15.79, 15.79, 15.79, 15.88, 15.88, 15.88, 15.88, 15.88, 15.91, 15.91, 15.91, 15.91, 15.91, 15.88, 15.88, 15.88, 15.88, 15.88, 15.91, 15.91, 15.91, 15.91, 15.91, 15.7, 15.7, 15.7, 15.7, 15.7, 15.79, 15.79, 15.79, 15.79, 15.79, 15.74, 15.74, 15.74, 15.74, 15.74, 15.81, 15.81, 15.81, 15.81, 15.81, 15.81, 15.81, 15.81, 15.81, 15.81, 15.87, 15.87 };
const int numData = sizeof(outputData) / sizeof(outputData[0]);

void printSOPDT() {
	StepFunction stepFunc(20.0, 30.0, 32.823277, 178, 1.0);
	SecondOrderPlusDelayLineOverdampedFunction func(&stepFunc, nullptr, 0);
	double params[4];
	params[PARAM_K] = 0.251778;
	params[PARAM_T] = 55.7078;
	params[PARAM_T2] = 55.7077;
	params[PARAM_L] = 1.80759;

	for (int i = 0; i < numData; i++) {
		double v = func.getEstimatedValueAtPoint(i, params);
		printf("%d,%f,%f,%f\r\n", i, outputData[i], v, stepFunc.getValue((float)i, 0));
	}
}

TEST(pidAutoTune, chsSopdtPid) {
	PidAutoTuneChrSopdt chr;
	double params0[4];
	
	// todo: find better initial values?
	params0[PARAM_K] = 0.1;
	params0[PARAM_T] = 1;
	params0[PARAM_T2] = 1;
	params0[PARAM_L] = 1;

	for (int tries = 0; tries < 10; tries ++) {
		for (int i = 0; i < numData; i++) {
			chr.addData(outputData[i]);
		}

		PidAutoTuneSettings settings;
		settings.minValue = 20.0;
		settings.maxValue = 30.0;
		settings.stepPoint = 178;	// todo: adjust for the buffer scale
		settings.maxPoint = 460;
		settings.timeScale = 1.0;
		if (chr.findPid(PID_TUNE_CHR2, settings, params0))
			break;
		// todo: the solver has failed. Choose other initial params?
		// params0[0] = ; params0[1] = ; params0[2] = ;
		break;
	}

#ifdef PID_DEBUG
	const double *p = chr.getParams();
	printf("Params: K=%g T1=%g T2=%g L=%g\r\n", p[PARAM_K], p[PARAM_T], p[PARAM_T2], p[PARAM_L]);
	const pid_s & pid = chr.getPid();
	printf("PID: P=%f I=%f D=%f offset=%f\r\n", pid.pFactor, pid.iFactor, pid.dFactor, pid.offset);
#endif

	// todo: check results
}

TEST(pidAutoTune, testPidCoefs) {
	pid_s pidParams[] = {
		{ 2.378598f, 0.011108f, 0.063678f, 32.823277f, 0, 100 },	// CHR1		ITAE=102.008 ISE=787.356 Overshoot=2.41016%
		{ 18.588152f, 0.166438f, 518.990417f, 32.823277f, 0, 100 },	// CHR2		ITAE=28.1444 ISE=186.409 Overshoot=7.86222%
		{ 2.383054f, 0.606178f, 0.067225f, 32.823277f, 0, 100 },	// CHR21	ITAE=215.15 Overshoot=18.5046
		{ 1.764889f, 0.106801f, 2.620852f, 32.823277f, 0, 100 },	// IMC21	ITAE=112.804 Overshoot=17.7643
		{ 292.501831f, 2.601279f, 8333.136719, 32.823277f, 0, 100 },// VDG2		ITAE=130.496 ISE=891.474 Overshoot=13.4626%
	};

	double params[4];
	params[PARAM_K] = 0.251778;
	params[PARAM_T] = 55.841;
	params[PARAM_T2] = 55.841;
	params[PARAM_L] = 1.52685;

	// todo: is it correct?
	double dTime = 1;
	for (int idx = 0; idx <= 4; idx++) {
		PidAccuracyMetric metric = PidAutoTuneChrSopdt::simulatePid<2048>(13.0, 14.0, dTime, pidParams[idx], params);
#ifdef PID_DEBUG
		printf("Metric result: ITAE=%g ISE=%g Overshoot=%g%%\r\n", metric.getItae(), metric.getIse(), metric.getMaxOvershoot() * 100.0);
#endif
	}

	// todo: check results
}

#if 0
GTEST_API_ int main(int argc, char **argv) {
	//testPidCoefs();

	//printSOPDT();
	testing::InitGoogleTest(&argc, argv);
	// uncomment if you only want to run selected tests
#if 0
	testMeasuredDataBuffer();
	testMatrixInverse();
	testGaussianFunction();
	testFOPDT();
	testSOPDTOverdamped();
	testChsSopdtPid();
#endif	
	::testing::GTEST_FLAG(filter) = "*testPidCoefs*";
	int result = RUN_ALL_TESTS();
	// windows ERRORLEVEL in Jenkins batch file seems to want negative value to detect failure
	return result == 0 ? 0 : -1;
}
#endif
