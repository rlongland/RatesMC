/* ======================================================================
   RatesMC2
   Author: R. Longland
   Start Date: 2018-11-13

   Description: Monte Carlo reaction rate code. This version 2 is for
   public release!

   TO DO:
    - Constants for resonances
    - Output S-factor for any broad resonances + non-resonant terms
    - Output broad resonance integrand? All resonances at all temperatures?
    - Output Porter-Thomas samples
    - Error accounting output
    - R codes for analysis
      - PlotUncertainty.R (like PlotCompare/PlotContour)
      - PlotPanel6.R
      - PlotPanelall.R
      - PlotSFactor.R (plots the S-factor for broad and non-resonant parts)
      - PlotIntegrand.R (plots the rate integrand at a given temperature)
      - PlotPT.R (plots the Porter-Thomas samples for a given resonance)
      - PlotCorrelations.R (plots some set of Rate vs. input parameter at a given T)
   ======================================================================
*/

#include <iostream>
#include <iomanip>
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

  // Input and output files
  std::string ofilename;
  std::string ofullfilename;
  std::string ifilename;
  if(argc==1){
    ifilename = "RatesMC.in";
    ofilename = "RatesMC.out";
    ofullfilename = "RatesMC.full";
  } else {
    std::cout << "Custom filenames not yet implemented!" << std::endl;
    return 0;
    ifilename = argv[1];
  }
  outfile.open(ofilename);
  outfullfile.open(ofullfilename);
  // Open log file
  logfile.open("RatesMC.log");
  // Open the test file for writing samples, etc.
  testfile.open("test.dat");
  // Open the file that stores Porter Thomas values
  ptfile.open("PT.dat");
  // LaTeX output file
  latexfile.open("RatesMC.latex");
  // Reaction rate sample
  sampfile.open("RatesMC.samp");
  // Contribution file
  contribfile.open("RatesMC.cont");

  // Make a reaction. This is where everything is held
  Reaction *Reac = new Reaction();
  
  // Open the input file
  int ret = ReadInputFile(ifilename, Reac);

  // Make output file headers
  writeOutputFileHeaders(Reac);
  Reac -> setupContribHeader();
  
  // Set up the random sampler
  setupRandom();
  
  // Prepare MC samples
  //  - For each resonance, sample all input parameters
  //  - Store every input parameter in a matrix (column = parameter, row = sample)
  //  - 
  Reac -> prepareSamples();

  // Write the reaction information to log file for diagnostics
  Reac -> writeReaction();

  // Write all samples to a file for later analysis
  Reac -> writeSamples();         

  // Loop through temperatures (this is parallelization happens)
  // At each temperature
  //  - Calculate rate
  //  - Store in rate matrix
  //  - 
  // First define the temperatures
  defineTemperatures();

  // Define the classical rates
  std::vector<double> classicalRate;

  // Now do the big loop over temperatures in parallel!!
  // Do all of the calculations first, then collect everything together
  omp_set_num_threads(1);
#pragma omp parallel for ordered
  for(double T : Temp){
    // ------------------------
    // FOR EACH TEMPERATURE
    // ------------------------
    
    int ID = omp_get_thread_num();
    std::cout << std::endl;
    std::cout << "--------------------------------------------------\n";
    std::cout << "Proc(" << ID << ") T = " << T; // << "\n";
    std::cout << "\n" ;

    logfile << "Temperature = " << T << " GK" << std::endl;
    
    // ------------------------
    // CALCULATE RATE
    // ------------------------
    
    // Calculate the non-resonant rate
    double ADRate[2];
    for(int j=0; j<2; j++){
      //      Old method using analytical formula
      //      ADRate[j] = Reac -> calcNonResonant(T, j);
      // New method of simply integrating the astrophysical s-factor
      ADRate[j] = Reac -> calcNonResonantIntegrated(T, j);
    }
    
    // Calculate the resonant rate
    double ResRate = Reac -> calcResonant(T);

    // ------------------------
    // COLLECT RATE
    // ------------------------

    // The classical rates can be easily summed
    classicalRate.push_back(ADRate[0]+ADRate[1]+ResRate);
    std::cout << "Classical Total Rate = " << classicalRate.back() << "\n";

    // Contribution array (NSamples)by(NRes+2)
    std::vector<std::vector<double> > Contributions;
    // Vector of rate samples at this temperature
    std::vector<double> RateSample;
    
    for(int s=0; s<NSamples; s++){

      // The non-resonant rate
      double ADRate0 = Reac -> getARate(s, 0);
      double ADRate1 = Reac -> getARate(s, 1);

      // Next get a vector that contains sample s from every resonance
      std::vector<double> resonancesSample = Reac -> getResonantRateSample(s);

      // Sum the total rate
      double totalRate = ADRate0 + ADRate1;
      for(double res : resonancesSample){
	totalRate += res;
      }

      // Calculate contribution for each resonance. 
      std::vector<double> Cont;
      Cont.push_back(ADRate0/totalRate);
      Cont.push_back(ADRate1/totalRate);
      for(double res : resonancesSample)
	Cont.push_back(res/totalRate);
      Contributions.push_back(Cont);
      
      // Fill the total reaction rate vector
      RateSample.push_back(totalRate);
      
    }

    // Write the contributions
    writeContributions(Contributions, T);

    // Write the rates + LaTeX
    writeRates(RateSample, classicalRate.back(), T);

    /*
#pragma omp critical
    {

      // print out the thread number and temperature
      // #pragma omp ordered
      std::cout << "Classical Resonant Rate (again) = " << ResRate << "\n";
    }
    */

    // Summarize any errors that may have occurred
    summarizeErrors(T);

  }

  
  
  std::cout << " ********************************************\n";
  std::cout << " *             Farewell!                    *\n";
  std::cout << " ********************************************\n";
  
  
  // Close the logfile
  logfile.close();
  testfile.close();
  ptfile.close();
  outfile.close();
  outfullfile.close();
  
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
