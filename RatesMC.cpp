/* RatesMC: Monte Carlo reaction rate code
 *
 * Copyright (C) 2022  Richard Longland
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: R. Longland
 *  Email: rllongla@ncsu.edu
 *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <cmath>

//#include "omp.h"

#include "gsl/gsl_sf_log.h"

#include "RatesMC.h"
#include "Utilities.h"
#include "Reaction.h"
#include "Resonance.h"


int main(int argc, char** argv){

  // Write the welcome screen
  WelcomeScreen();

  //  gsl_error_handler_t myHandler;
  //  gsl_set_error_handler(&);
  //gsl_sf_log(0.0);
  //gsl_set_error_handler(NULL);
  
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
  integrandfile.open("RatesMC.integ");
  // Open the file that stores Porter Thomas values
  ptfile.open("RatesMC.PT");
  // LaTeX output file
  latexfile.open("RatesMC.latex");
  // Reaction rate sample
  sampfile.open("RatesMC.samp");
  // Contribution file
  contribfile.open("RatesMC.cont");
	// Test file for storing debugging
	testfile.open("test.dat");
	
  // Make a reaction. This is where everything is held
  Reaction *Reac = new Reaction();
  
  // Open the input file
  int ret = ReadInputFile(ifilename, Reac);

	if(ret != 0){
		std::cout << "ERROR: You should never see this message!\n";
		exit(EXIT_FAILURE);
	}
	
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

	// Before we write anything long, write the astrophysical S-factor
	Reac -> writeSFactor();

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
	std::vector<double>::iterator it;
  //omp_set_num_threads(1);
	//#pragma omp parallel for ordered
	for(it = Temp.begin(); it < Temp.end(); ++it){
    // ------------------------
    // FOR EACH TEMPERATURE
    // ------------------------
		double T = *it;
    
    int ID = 0;//omp_get_thread_num();
    std::cout << std::endl;
    std::cout << "--------------------------------------------------\n";
    std::cout << "Proc(" << ID << ") T = " << T; // << "\n";
    std::cout << "\n" ;
    

    
    //logfile << "--------------------\n";
    logfile << "Temperature = " << T << " GK" << std::endl;
		integrandfile << "Temperature = " << T << " GK" << std::endl;
		
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
    //std::cout << "Classical Total Rate = " << classicalRate.back() << "\n";

		// Combine Resonance possibilities into the main corresponding resonance
		Reac -> CombineResonancePossibilities();
		
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
				//if( !std::isnan(res) )
				totalRate += res;
      }
			//std::cout << s << "  " << totalRate << "\n";

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

    // Write the rate samples
    writeRateSamples(RateSample, T);
    
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
		integrandfile << std::endl;
  }

  
  
  std::cout << " **************************************************\n";
  std::cout << " *                 Farewell!                      *\n";
  std::cout << " **************************************************\n";
  
  
  // Close the logfile
  logfile.close();
  integrandfile.close();
  ptfile.close();
  outfile.close();
  outfullfile.close();
	testfile.close();
	
  return 1;
}



void WelcomeScreen(){
  std::cout << std::endl;
  std::cout << " **************************************************" << std::endl;
  std::cout << " *                   RatesMC                      *" << std::endl;
	std::cout << " *        Copyright (C) 2022  R. Longland         *" << std::endl;
  std::cout << " *            V. " << VersionNumber << "  " << VersionDate 
						<< "                *" << std::endl;
	std::cout << " * This program comes with ABSOLUTELY NO WARRANTY *" << std::endl;
  std::cout << " *                                                *" << std::endl;
  std::cout << " **************************************************" << std::endl;
  std::cout << "\n" << std::endl;

  return;
}
