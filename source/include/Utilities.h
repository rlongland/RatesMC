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
  NANCount,BelowIntLimit,IntfNANCount;

// Histogram and output
extern double HistT, HistMin, HistMax;
extern std::ofstream logfile;
extern std::ofstream testfile;

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
