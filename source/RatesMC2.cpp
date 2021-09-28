/* ======================================================================
   RatesMC2
   Author: R. Longland
   Date: 2018-11-13

   Description: Monte Carlo reaction rate code. This version 2 is for
   public release!
   ======================================================================
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <algorithm>    // std::sort

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

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
  
  // Set up the random sampler
  setupRandom();
  
  // Prepare MC samples
  //  - For each resonance, sample all input parameters
  //  - Store every input parameter in a matrix (column = parameter, row = sample)
  //  - 
  Reac -> prepareSamples();

  // Write the reaction information to log file for diagnostics
  Reac -> writeReaction();
  //Reac -> printReaction();
  
  Reac -> writeSamples();         // Write all samples to a file

  // Loop through temperatures (this is parallelization happens)
  // At each temperature
  //  - Calculate rate
  //  - Store in rate matrix
  //  - 
  // First define the temperatures
  defineTemperatures();

  // Define the classical rates
  std::vector<double> classicalRate;

  // Open the output files
  sampfile.open("RatesMC.samp");
  contribfile.open("RatesMC.cont");

  Reac -> setupContribHeader();
  
  // Now do the big loop over temperatures in parallel!!
  // Do all of the calculations first, then collect everything together
  omp_set_num_threads(1);
#pragma omp parallel for ordered
  for(double T : Temp){
    int ID = omp_get_thread_num();
    std::cout << std::endl;
    std::cout << "--------------------------------------------------\n";
    std::cout << "Proc(" << ID << ") T = " << T; // << "\n";
    std::cout << "\n" ;

    // --------------
    // CALCULATE RATE
    // --------------
    // Calculate the non-resonant rate
    double ADRate[2];
    for(int j=0; j<2; j++){
      ADRate[j] = Reac -> calcNonResonant(T, j);
      std::cout << "ADRate = " << ADRate[j] << "\n";
    }
    
    // Calculate the resonant rate
    double ResRate = Reac -> calcResonant(T);

    // --------------
    // COLLECT RATE
    // --------------
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
	std::cout << res << " ";
	totalRate += res;
      }
      std::cout << "\n";

      std::cout << "Total rate = " << totalRate << "\n";

      // Calculate contribution for each resonance. 
      std::vector<double> Cont;
      Cont.push_back(ADRate0/totalRate);
      Cont.push_back(ADRate1/totalRate);
      for(double res : resonancesSample)
	Cont.push_back(res/totalRate);
      for(int i=0; i<Cont.size(); i++)
	std::cout << Cont[i] << " ";
      std::cout << "\n";
      Contributions.push_back(Cont);
      
      // Fill the total reaction rate vector
      RateSample.push_back(totalRate);
      
    }
    
    std::cout << "At the end we have\n";
    std::cout << "RateSample of length: " << RateSample.size() << "\n";
    //transpose(Contributions);
    std::cout << "Contributions of length: " << Contributions.size() << " X "
	      << Contributions[0].size() << "\n";
    std::cout << std::endl;
    
    // From the contributions, calculate Low, median, and high
    // contribution for each resonance
    // Sort, and find quantiles for each rate contribution
    std::vector<double> LowCont;
    std::vector<double> HighCont;
    std::vector<double> MedianCont;
    for(int j=0;j<Contributions.size();j++){

      // Sort the Contributions for each resonance
      //std::sort (Contributions[j].begin(), Contributions[j].end());
      for(int s = 0; s<Contributions[j].size(); s++){
	std::cout << Contributions[j][s] << " ";
      }
      std::cout << std::endl;
      
      // Now convert this into a gsl vector
      gsl_vector_const_view gsl_v =
	gsl_vector_const_view_array( &Contributions[j][0], Contributions[j].size() );
      
      LowCont.push_back(gsl_stats_quantile_from_sorted_data(gsl_v.vector.data,1, Contributions[j].size(),0.16));
      HighCont.push_back(gsl_stats_quantile_from_sorted_data(gsl_v.vector.data,1, Contributions[j].size(),0.84));
      MedianCont.push_back(gsl_stats_median_from_sorted_data(gsl_v.vector.data,1,Contributions[j].size()));
    }

    for(int i=0; i<LowCont.size(); i++){
      std::cout << "i = " << i << "  LowCont = " << LowCont[i] <<  "  MedianCont = " << MedianCont[i]
		<< "  HighCont = " << HighCont[i] << "\n";
    }
    
    /*
#pragma omp critical
    {

      // print out the thread number and temperature
      // #pragma omp ordered
      std::cout << "Classical Resonant Rate (again) = " << ResRate << "\n";
    }
    */
    // The classical rates can be easily summed
    classicalRate.push_back(ADRate[0]+ADRate[1]+ResRate);
    std::cout << "Classical Total Rate = " << classicalRate.back() << "\n";
  }

  
  
  std::cout << " ********************************************\n";
  std::cout << " *             Farewell!                    *\n";
  std::cout << " ********************************************\n";
  
  
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
