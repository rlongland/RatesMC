#ifndef _Utilities_h_
#define _Utilities_h_

#include "Reaction.h"
#include <gsl/gsl_rng.h>

// Random number generator
extern gsl_rng * r;

// program control parameters
extern int NSamples,NHists;
extern int NTemps;
extern double EMin,R0;
extern double EPS;

// Include correlations
extern bool bPartialWidthCorrelations;
extern bool bEnergyCorrelations;

// Counters
extern int PenZeroCount,IntegratedCount,SubSampledPosCount,SampledNegCount,
  NANCount,BelowIntLimit,IntfNANCount, LogZeroCount;

// Flags
extern bool ErrorFlag;

// Histogram and output
extern double HistT, HistMin, HistMax;
extern std::ofstream logfile;
extern std::ofstream testfile;
extern std::ofstream ptfile;
extern std::ofstream sampfile;
extern std::ofstream contribfile;
extern std::ofstream outfile;
extern std::ofstream outfullfile;
extern std::ofstream latexfile;

// Other program-wide variables
extern std::vector<double> Temp;

// Read input file
int ReadInputFile(std::string inputfilename, Reaction *R);

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

// Check how well the lognormal fits
double CalcAD(std::vector<double> Rates,double Mu,double Sigma);

// Check for zero
bool isZero(double x);

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
