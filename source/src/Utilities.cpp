#include <iostream>
#include <fstream>
#include <limits>

#include "Utilities.h"

double EPS=1.0e-5;


// Read an integer from a single line
int readInt(std::ifstream &infile){

  
  int input;
  infile >> input;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  return input;

}
// Read a double from a single line
double readDouble(std::ifstream &infile){

  double input;
  infile >> input;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  return input;

}
void skipLines(std::ifstream &infile, int nlines){

  std::string dummy;
  for(int i=0; i<nlines; i++){
    infile >> dummy;
    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
  
}

int ReadInputFile(std::string inputfilename, Reaction *R){

  std::cout << "The input file name is: " << inputfilename << std::endl;

  // Start by opening input file (reactions.dat)
  std::ifstream infile;
  infile.open(inputfilename.c_str(),std::ifstream::in);

  // test to see if file is there
  if(!infile.is_open()){
    std::cout << "ERROR! File: RatesMC.in is not found!" << std::endl;
    return(1);
  }

  std::string name,dummy;
  infile >> name;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setName(name);
  R -> getName();

  // Blank line
  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  // Masses, charges, spins, Q-values at the top
  int z0,z1,z2;
  double R0;
  double m0,m1,m2;
  double j0,j1,j2;
  double Qin,Qout;
  int gindex;
  
  infile >> z0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> z1;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> z2;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setCharges(z0,z1,z2);

  infile >> m0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> m1;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> m2;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setMasses(m0,m1,m2);

  infile >> j0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> j1;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> j2;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setSpins(j0,j1,j2);

  // Entrance and exit particle separation energies
  Qin = readDouble(infile);
  Qout = readDouble(infile);
  R -> setSeparationEnergies(Qin, Qout);

  // The R0 radius
  R0 = readDouble(infile);
  R -> setR0(R0);

  // Index of the gamma-ray channel
  gindex = readInt(infile);
  R -> setGammaIndex(gindex);

  // Ignore a line
  skipLines(infile, 1);

  // Computation control block
  int EMin = readDouble(infile);
  int NSamples = readInt(infile);
  int NTemps = readInt(infile);
  int itmp;
  infile >> itmp;
  bool bPartialWidthCorrelations = (bool)itmp;
  infile >> itmp;
  bool bEnergyCorrelations = (bool)itmp;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  
  return 0;
}

