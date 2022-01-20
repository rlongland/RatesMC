/* RatesMC: Monte Carlo reaction rate code
 *
 * Copyright (C) 2022  Richard Longland
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: R. Longland
 *  Email: rllongla@ncsu.edu
 *
 */

#ifndef _Utilities_h_
#define _Utilities_h_

#include "Reaction.h"
#include <gsl/gsl_rng.h>

extern std::string VersionNumber;
extern std::string VersionDate;

// Random number generator
extern gsl_rng * r;

// program control parameters
extern int NSamples,NHists;
extern int NTemps;
extern double EMin,R0;
extern double EPS;

extern double ElectronMass;
extern double AMU;

// Include correlations
extern bool bPartialWidthCorrelations;
extern bool bEnergyCorrelations;

// Counters
extern int PenZeroCount,IntegratedCount,SubSampledPosCount,SampledNegCount,
  NANCount,InfCount,BelowIntLimit,IntfNANCount, LogZeroCount;

// Flags
extern bool ErrorFlag;

// Histogram and output
extern double HistT, HistMin, HistMax;
extern std::ofstream logfile;
extern std::ofstream integrandfile;
extern std::ofstream ptfile;
extern std::ofstream sampfile;
extern std::ofstream contribfile;
extern std::ofstream outfile;
extern std::ofstream outfullfile;
extern std::ofstream latexfile;
extern std::ofstream testfile;

// Other program-wide variables
extern std::vector<double> Temp;

// Read input file
int ReadInputFile(std::string inputfilename, Reaction *R);

// Is a read thing a number?
bool isNumber(const std::string& str);
// Read integer and double from single line
int readInt(std::ifstream &in);
double readDouble(std::ifstream &in);
// Skip n lines
void skipLines(std::ifstream &in, int n);
// Count entries in a line
int countEntries(std::ifstream &in);
// Read a non-resonant line
void readNonResonant(std::ifstream &in, Reaction &R, int part);
// Read all of the standard resonances
void readResonanceBlock(std::ifstream &in, Reaction &R);

// Define the temperature array
void defineTemperatures();

// Utility to take expectation value and variance and turn them into lognormal mu and sigma
void logNormalize(double exp, double sd, double& mu, double& sigma);

// Penetration Factor calculation
double PenFactor(double E, double L, double Mass0, double Mass1,
		 int Charge0, int Charge1, double R);

// Transpose a 2D double vector
void transpose(std::vector<std::vector<double> > &b);

// Write headers (reaction name, etc) into output files
void writeOutputFileHeaders(Reaction *R);
// Write the contributions of each resonance to file
void writeContributions(std::vector<std::vector<double> > Contributions, double Temperature);
// Write the rates to file
void writeRates(std::vector<double> Rates, double ARate, double Temperature);
// Write the LaTeX file
void WriteLatex2(double Temperature, double LowRate, double MedianRate, double HighRate,
		 double RateSigma);
// Write the sample file
void writeRateSamples(std::vector<double> RateSample, double Temp);
// Write a summary of all errors at this temperature to the log file
void summarizeErrors(double Temp);

// Check how well the lognormal fits
double CalcAD(std::vector<double> Rates,double Mu,double Sigma);

// Check for zero
bool isZero(double x);

// Convert atomic to nuclear mass
double atomicToNuclear(double A, double Z);

// Check whether a file exists
inline bool fileExists (const std::string& name);

// Setup the random number generator
void setupRandom();
unsigned long int random_seed();

// Count the number of numbers in a single line
class NumberCounter{
 public:
  NumberCounter();
  void count(std::ifstream &in);
  friend std::ostream& operator<<(std::ostream &out, const NumberCounter &w);

 private:
  int numbers;
  int spaces;
  int comments;
};


#endif
