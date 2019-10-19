/*
* @file	pid_from_msl.cpp
*
* @date Oct 5, 2019
* @author andreika, (c) 2019
*/

#include "global.h"
#include <fstream>
#include <vector>

#include "pid_auto.h"


class MslData {
public:
	bool readMsl(const char *fname, double startTime, double endTime, int inputIdx, int outputIdx) {
		std::ifstream fp(fname);

		if (!fp)
			return false;
		
		curIdx = -1;
		settings.minValue = settings.maxValue = settings.maxPoint = 0;
		settings.timeScale = 1,
		settings.stepPoint = -1.0;
		totalTime = 0;

		std::string str;
		for (int i = 0; std::getline(fp, str); i++) {
			// data starts on the 4th line
			if (i < 4)
				continue;
			parseLine(str, startTime, endTime, inputIdx, outputIdx);
		}
		
		settings.maxPoint = getSaturationStartPoint();
		assert(data.size() == curIdx);

		settings.timeScale = settings.maxPoint / totalTime;

		fp.close();
		return true;
	}

	bool parseLine(const std::string & str, double startTime, double endTime, int inputIdx, int outputIdx) {
		std::stringstream sstr(str);
		std::string item;
		for (int j = 0; getline(sstr, item, '\t'); j++) {
			double v = atof(item.c_str());
			// the first column is timestamp
			if (j == 0) {
				if (v < startTime || v > endTime)
					return false;
				if (curIdx < 0)
					prevTime = v;
				totalTime += v - prevTime;
				prevTime = v;
			} else if (j == inputIdx) {
				// this is an input step, we should find it
				if (curIdx < 0) {
					settings.minValue = v;
					curIdx = 0;
				} else if (v != settings.minValue && settings.stepPoint < 0) {
					settings.maxValue = v;
					settings.stepPoint = curIdx;
				}
				curIdx++;
			} else if (j == outputIdx) {
				data.push_back((float)v);
				if (curIdx >= 0 && settings.stepPoint < 0) {
					// calculate averaged level to determine the acceptable noise level
					averagedMin = (averagedMin * (curIdx - 1) + v) / curIdx;
					// this is not accurate because 'averagedMin' is continuously changing
					acceptableNoiseLevel = std::max(acceptableNoiseLevel, abs(v - averagedMin));
				}
			}
		}
		return true;
	}

	double getSaturationStartPoint() {
		int i;
		double j;
		// max noise level is used to get the saturation limit of the signal
		double curNoiseLevel = 0, averagedMax = 0;
		// we step back some points from the last one and find the saturation start
		for (i = curIdx - 1, j = 1.0; i > settings.stepPoint; i--, j += 1.0) {
			double v = data[i];
			averagedMax = (averagedMax * (j - 1) + v) / j;
			// this is not accurate because 'averagedMax' is continuously changing
			curNoiseLevel = std::max(curNoiseLevel, abs(v - averagedMax));
			
			// we assume that the "upper" level noise is like the same as the "lower" noise, so we compare them,
			// and if the noise level starts growing, then we're in the step transient zone, and we stop
			if (curNoiseLevel > acceptableNoiseLevel) {
				break;
			}
		}
		// return the point in the middle, just to be safe (we don't want to be close to the transient zone)
		return (curIdx - 1 + i) / 2;
	}

public:
	std::vector<float> data;
	double totalTime = 0, prevTime = 0;
	PidAutoTuneSettings settings;
	int curIdx = -1;
	float prevV = 0;
	
	// we assume that the signal is quasi-stable (asymptotic) from the start point until the 'stepPoint';
	// it's noise level is used to find the saturation limit of the rest of the data (see getSaturationStartPoint())
	double acceptableNoiseLevel = 0;
	double averagedMin = 0;
};

const char *getMethodName(pid_tune_method_e method) {
	switch (method) {
	case PID_TUNE_CHR1: return "CHR1";
	case PID_TUNE_AUTO1: return "AUTO1";
	case PID_TUNE_IMC2_1: return "IMC2_1";
	case PID_TUNE_CHR2_1: return "CHR2_1";
	case PID_TUNE_CHR2: return "CHR2";
	case PID_TUNE_VDG2: return "VDG2";
	case PID_TUNE_HP2: return "HP2";
	case PID_TUNE_AUTO2: return "AUTO2";
	}
	return "(N/A)";
}

const char *getSimTypeName(pid_sim_type_e simType) {
	switch (simType) {
	case PID_SIM_REGULATOR: return "Regulator";
	case PID_SIM_SERVO: return "Servo";
	}
	return "(N/A)";
}

#if 1
int main(int argc, char **argv) {
	if (argc < 6) {
		printf("Usage: PID_FROM_MSL file.msl start_time end_time input_column output_column...\r\n");
		return -1;
	}
	
	printf("PID_FROM_MSL - find PID controller coefficients based on a measured step response in a rusEFI log file.\r\n");
	printf("Version 0.2 (c) andreika, 2019\r\n\r\n");
	printf("Reading file %s...\r\n", argv[1]);
	
	MslData data;
	if (!data.readMsl(argv[1], atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]))) {
		printf("Usage: PID_FROM_MSL <file.msl> <startTime> <endTime> <inColumnIdx> <outColumnIdx> [<targetValue>]\r\n");
		return -2;
	}

	// Target value is optional, for PID_SIM_REGULATOR only
	data.settings.targetValue = (argc > 6) ? atof(argv[6]) : 0.0;

	printf("Measuring Settings: targetValue=%g minValue=%g maxValue=%g stepPoint=%g maxPoint=%g numPoints=%d timeScale=%g\r\n",
		data.settings.targetValue, data.settings.minValue, data.settings.maxValue,
		data.settings.stepPoint, data.settings.maxPoint, data.data.size(), data.settings.timeScale);

	PidAutoTune chr1, chr2;
	static const int numPids = 2;
	PidAutoTune *chr[numPids] = { &chr1, &chr2 };
	
	for (size_t i = 0; i < data.data.size(); i++) {
		for (int j = 0; j < numPids; j++) {
			chr[j]->addData(data.data[i]);
		}
	}

	// todo: more flexible method chooser
	pid_sim_type_e simTypes[numPids] = { PID_SIM_REGULATOR, PID_SIM_REGULATOR };
	pid_tune_method_e methods[numPids] = { PID_TUNE_AUTO1, PID_TUNE_AUTO2 };
	
	for (int k = 0; k < numPids; k++) {
		printf("\r\n%d) Trying method %s on \"%s\" PID model:\r\n", k + 1, getMethodName(methods[k]), getSimTypeName(simTypes[k]));
		chr[k]->findPid(simTypes[k], methods[k], data.settings, nullptr);
	}

	printf("Done!\r\n");

	// todo: is it correct?
	double dTime = 1.0 / data.settings.timeScale;
	const int numSimPoints = 1024;

	pid_s bestPid;
	double smallestItae = DBL_MAX;
	for (int k = 0; k < 2; k++) {
		const double *p = chr[k]->getParams();
		printf("Model-%d Params: K=%g T1=%g T2=%g L=%g\r\n", (k + 1), p[PARAM_K], p[PARAM_T], p[PARAM_T2], p[PARAM_L]);

		pid_s pid[2] = { chr[k]->getPid0(), chr[k]->getPid() };
		for (int j = 0; j < 2; j++) {

			const char *csvName = (j == 0) ? ((k == 0) ? "pid_test01.csv" : "pid_test02.csv") : ((k == 0) ? "pid_test1.csv" : "pid_test2.csv");
			
			PidSimulator<numSimPoints> sim1(simTypes[k], chr[k]->getMethodOrder(methods[k]),
					chr[k]->getAvgMeasuredMin(), chr[k]->getAvgMeasuredMax(), dTime, chr[k]->getModelBias(), csvName);

			PidAccuracyMetric metric = sim1.simulate(numSimPoints, pid[j], p);
			
			printf("  PID%d: P=%.8f I=%.8f D=%.8f offset=%.8f period=%.8fms\r\n", j, pid[j].pFactor, pid[j].iFactor, pid[j].dFactor, pid[j].offset,
				pid[j].periodMs);
			printf("  Metric%d result: ITAE=%g ISE=%g Overshoot=%g%%\r\n", j, metric.getItae(), metric.getIse(), metric.getMaxOvershoot() * 100.0);

			if (metric.getItae() < smallestItae) {
				smallestItae = metric.getItae();
				bestPid = pid[j];
			}
		}
	}

	printf("The best PID: P=%.8f I=%.8f D=%.8f offset=%.1f period=%.1fms\r\n", bestPid.pFactor, bestPid.iFactor, bestPid.dFactor, bestPid.offset, bestPid.periodMs);

	return 0;
}
#endif
