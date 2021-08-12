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

#include "omp.h"

#include "RatesMC2.h"
#include "Utilities.h"
#include "Reaction.h"
#include "Resonance.h"
#include "DirectCapture.h"


int main(int argc, char** argv){

  // Write the welcome screen
  WelcomeScreen();

  // Open log file
  logfile.open("RatesMC.log");

  // Open the test file for writing samples, etc.
  testfile.open("test.dat");

  // Open the file that stores Porter Thomas values
  ptfile.open("PT.dat");
  
  // Input and output files
  std::string ofilename;
  std::string ifilename;
  if(argc==1){
    ifilename = "RatesMC.in";
    ofilename = "RatesMC.out";
  } else {
    ifilename = argv[1];
  }
  
  // Make a reaction. This is where everything is held
  Reaction *Reac = new Reaction();
  
  // Open the input file
  int ret = ReadInputFile(ifilename, Reac);
  // Write the reaction information to log file for diagnostics
  Reac -> writeReaction();

  // Set up the random sampler
  setupRandom();
  
  // Prepare MC samples
  //  - For each resonance, sample all input parameters
  //  - Store every input parameter in a matrix (column = parameter, row = sample)
  //  - 
  Reac -> prepareSamples();


  // Loop through temperatures (this is parallelization happens)
  // At each temperature
  //  - Calculate rate
  //  - Store in rate matrix
  //  - 
  // First define the temperatures
  defineTemperatures();
  logfile << "--------------------------------------------------\n";
  logfile << "There are " << Temp.size() << " temperatures:\n";

  // Now do the big loop in parallel!!
  //omp_set_num_threads(4);
#pragma omp parallel for ordered
  for(double T : Temp){
    int ID = omp_get_thread_num();

    // print out the thread number and temperature
    // #pragma omp ordered
#pragma omp critical
    //std::cout << "(" << ID << ") " << T << "\n";
    std::cout << " " ;
  }
  
  std::cout << "\n";

  
  
  // Close the logfile
  logfile.close();
  testfile.close();
  ptfile.close();
  
  return 1;
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
