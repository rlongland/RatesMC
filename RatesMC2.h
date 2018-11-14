#ifndef _RatesMC2_h_
#define _RatesMC2_h_

#include <fstream>
#include <iostream>

#include "Reaction.h"

using std::string;
using std::cout;
using std::endl;

// The welcome screen
void WelcomeScreen();
// Read integer and double from single line
int readInt();
double readDouble();
// Read input file
int ReadInputFile(std::string inputfilename, Reaction *R);


// Versions
std::string VersionNumber = "2.0.0";
std::string VersionDate = "Nov. 13, 2018  ";

// program control parameters
int NSamples,NHists;
int NTemps;
double CutoffE[2],EMin,R0;
double EPS=1.0e-5;

// Include correlations
bool bPartialWidthCorrelations;
bool bEnergyCorrelations;

// Counters
int PenZeroCount,IntegratedCount,SubSampledPosCount,SampledNegCount,
  NANCount,BelowIntLimit,IntfNANCount;

// Histogram and output
double HistT, HistMin, HistMax;
std::ofstream logfile;
std::ofstream testfile;

#endif
