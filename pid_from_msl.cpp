/*
* @file	pid_from_msl.cpp
*
* @date Oct 5, 2019
* @author andreika, (c) 2019
*/

#include "global.h"
#include <fstream>
#include <vector>

#include "pid_open_loop_models.h"
#include "pid_controller.h"


class MslData {
public:
	bool readMsl(const char *fname, double startTime, double endTime, int inputIdx, int outputIdx) {
		std::ifstream fp(fname);

		if (!fp)
			return false;
		
		curIdx = -1;
		stepPoint = -1.0;

		std::string str;
		for (int i = 0; std::getline(fp, str); i++) {
			// data starts at 4th line
			if (i < 4)
				continue;
			parseLine(str, startTime, endTime, inputIdx, outputIdx);
		}
		
		maxPoint = curIdx - 1;
		assert(data.size() == curIdx);

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
			} else if (j == inputIdx) {
				// this is an input step, we should find it
				if (curIdx < 0) {
					minValue = v;
					curIdx = 0;
				} else if (v != minValue && stepPoint < 0) {
					maxValue = v;
					stepPoint = curIdx;
				}
				curIdx++;
			} else if (j == outputIdx) {
				data.push_back((float)v);
			}
		}
		return true;
	}

public:
	std::vector<float> data;
	double minValue = 0, maxValue = 0, stepPoint = -1.0, maxPoint = 0;
	int curIdx = -1;
};

#if 0
int main(int argc, char **argv) {
	if (argc < 6) {
		printf("Usage: PID_FROM_MSL file.msl start_time end_time input_column output_column...\r\n");
		return -1;
	}
	
	printf("PID_FROM_MSL: Reading file %s...\r\n", argv[1]);
	
	MslData data;
	if (!data.readMsl(argv[1], atof(argv[2]), atof(argv[3]), atoi(argv[4]), atoi(argv[5]))) {
		return -2;
	}

	printf("PID_FROM_MSL: Calculating...\r\n");

	PidAutoTuneChrSopdt chr;
	double params0[4];

	// todo: find better initial values?
	params0[PARAM_K] = 0.1;
	params0[PARAM_T] = 1;
	params0[PARAM_T2] = 1;
	params0[PARAM_L] = 1;

	for (size_t i = 0; i < data.data.size(); i++) {
		chr.addData(data.data[i]);
	}

	double stepPoint = 178;	// todo: adjust for the buffer scale
	double maxPoint = 460;
	chr.findPid(PID_TUNE_CHR2, data.minValue, data.maxValue, data.stepPoint, data.maxPoint, params0);

	printf("Done!\r\n");

	const double *p = chr.getParams();
	printf("Model Params: K=%g T1=%g T2=%g L=%g\r\n", p[PARAM_K], p[PARAM_T], p[PARAM_T2], p[PARAM_L]);
	const pid_s & pid = chr.getPid();
	printf("PID: P=%f I=%f D=%f offset=%f\r\n", pid.pFactor, pid.iFactor, pid.dFactor, pid.offset);


	return 0;
}
#endif
