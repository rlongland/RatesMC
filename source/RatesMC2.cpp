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
#include "Utilities.h"
#include "Reaction.h"
#include "Resonance.h"
#include "DirectCapture.h"


int main(int argc, char** argv){

  // Write the welcome screen
  WelcomeScreen();

  // Log file
  //std::ofstream logfile;
  logfile.open("RatesMC.log");
  
  // Input and output files
  std::string ofilename;
  std::string ifilename;
  if(argc==1){
    ifilename = "RatesMC.in";
    ofilename = "RatesMC.out";
  } else {
    ifilename = argv[1];
  }

  // Keep all of the settings in a settings class
  // Settings Set;
  
  // Make a reaction. This is where everything is held
  Reaction *Reac = new Reaction();
  
  // Open the input file
  int ret = ReadInputFile(ifilename, Reac);
  
  Reac -> writeReaction();
  
  /*
    Testing of a few things
  */
  DirectCapture DC;

  logfile.close();
  
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
