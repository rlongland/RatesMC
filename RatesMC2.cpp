/* ======================================================================
   RatesMC2
   Author: R. Longland
   Date: 2018-11-13

   Description: Monte Carlo reaction rate code. This version 2 is for
   public release!
   ======================================================================
*/

#include <iostream>
#include <fstream>
#include <limits>

#include "RatesMC2.h"
#include "Reaction.h"
#include "Resonance.h"
#include "DirectCapture.h"


int main(int argc, char** argv){

  // Write the welcome screen
  WelcomeScreen();

  string ofilename;
  string ifilename;
  
  if(argc==1){
    ifilename = "RatesMC.in";
    ofilename = "RatesMC.out";
  } else {
    ifilename = argv[1];
  }

  // Make a reaction. This is where everything is held
  Reaction Reac;
  
  // Open the input file
  int ret = ReadInputFile(ifilename, &Reac);
  
  Reac.printReaction();
  
  /*
    Testing of a few things
  */
  DirectCapture DC;

  return 1;
}

int ReadInputFile(std::string inputfilename, Reaction *R){

  cout << "The input file name is: " << inputfilename << endl;

  // Start by opening input file (reactions.dat)
  std::ifstream infile;
  infile.open(inputfilename.c_str(),std::ifstream::in);

  // test to see if file is there
  if(!infile.is_open()){
    cout << "ERROR! File: RatesMC.in is not found!" << endl;
    return(1);
  }

  string name,dummy;
  infile >> name;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setName(name);
  R -> getName();

  // Blank line
  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  // Masses, charges, spins, Q-values at the top
  int z0,z1,z2;
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

  infile >> Qin;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> Qout;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setSeparationEnergies(Qin, Qout);

  infile >> R0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  infile >> gindex;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setGammaIndex(gindex);

  // Ignore a line
  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  
  
  return 1;
}

// Read an integer from a single line
int readInt(){

  /*
  string input;
  infile >> input;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  int ret = (int)input;
  */
  int ret = 0;
  return ret;

}
// Read a double from a single line
double readDouble(){

  /*
  string input;
  infile >> input;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  double ret = (double)input;
  cout << "Read " << ret << endl;
  */
  double ret = 0;
  return ret;

}

void WelcomeScreen(){
  std::cout << std::endl;
  std::cout << " ********************************************" << std::endl;
  std::cout << " *           Welcome to RatesMC             *" << std::endl;
  std::cout << " *         V. " << VersionNumber << "  " << VersionDate 
       << "        *" << std::endl;
  std::cout << " *     Written by:   Richard Longland       *" << std::endl;
  std::cout << " *                                          *" << std::endl;
  std::cout << " ********************************************" << std::endl;
  std::cout << "\n" << std::endl;

  return;
}
