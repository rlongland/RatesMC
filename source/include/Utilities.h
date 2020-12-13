#ifndef _Utilities_h_
#define _Utilities_h_

#include "Reaction.h"

  
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
