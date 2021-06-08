#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <sstream>

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

// Count the number of entries on a line
int countEntries(std::ifstream &infile){

  std::string line;
  std::vector<std::string> entries;

  // First save position
  int place=infile.tellg();

  std::getline(infile, line);
  std::stringstream ss(line);
  std::string entry;
  while( ss >> entry ){
    entries.push_back(entry);
    //std::cout << entry << " ";
  }
  //std::cout << "\n";

  // Return to saved position in file
  infile.seekg(place);
  
  return entries.size();
  
}

void readNonResonant(std::ifstream &infile, Reaction &R, int part){

  int ne = countEntries(infile);
  std::cout << "There are " << ne << " entries\n";
  if(ne != 5){
    std::cout << "  ERROR: There should be 5 numbers in the non-resonant input lines\n";
    exit(EXIT_FAILURE);
  }

  double s, sp, spp, ds, cutoffe;
  infile >> s >> sp >> spp >> ds >> cutoffe;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R.setNonResonant(s,sp,spp,ds,cutoffe,part);

}

/* 
   This is where Victor need to work on reading all of the "standard" resonances
   1) Start by just reading 5 of them
   2) Then add the ability to add up until the end of the block (check for ****?)
   3) Include error checking. Do the numbers make sense?
*/

void readResonanceBlock(std::ifstream &infile, Reaction &R){

  double E_cm, dE_cm, wg, dwg;
  
  for(int i=0; i<5; i++){

    infile >> E_cm >> dE_cm >> wg >> dwg;
    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    std::cout << E_cm << " " << dE_cm << " " << wg << " " << dwg << "\n";
    R.addResonance(i, E_cm, dE_cm, wg, dwg);

  }

}

int ReadInputFile(std::string inputfilename, Reaction *R){

  std::cout << "The input file name is: " << inputfilename << std::endl;

  int ne;  // for counting entries
  
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

  // Skip 3 lines
  skipLines(infile, 3);

  // Non-resonant line 1
  readNonResonant(infile, *R, 0);
  
  // Non-resonant line 2
  readNonResonant(infile, *R, 1);

  // Skip 5 lines
  skipLines(infile, 5);

  // Read all of the resonances <- VICTOR'S PROJECT
  readResonanceBlock(infile, *R);
  
  return 0;
}

