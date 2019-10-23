/*
* @file	output_csv.h
*
* Output data to .CSV file
*
* @date Oct 06, 2019
* @author andreika, (c) 2019
*/

#pragma once

// Use variadic templates for output

static std::ofstream csvfp;
static const char *csvfname = nullptr;

template<typename T>
static void output_csv(const char *fname, T arg) {
	if (fname == nullptr)
		csvfp << ";";
	csvfp << arg;
	//if (fname != nullptr)
	//	csvfp << std::endl;
}

template<typename T, typename... Args>
static void output_csv(const char *fname, T v, Args... args) {
	if (fname != csvfname && fname != nullptr) {
		csvfp.close();
		csvfp.open(fname);
		csvfname = fname;
	}
	output_csv(fname, v);
	output_csv(nullptr, args...);
	if (fname != nullptr)
		csvfp << std::endl;
}


