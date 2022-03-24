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
#include <vector>
#include <sstream>
#include <algorithm>    // std::sort
#include <cstdlib>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

#include "Utilities.h"
#include "AMEReader.h"

double EPS=1.0e-5;
double ElectronMass=0.0005485799;   // in AMU
double AMU = 931.494102;            // AMU -> MeV
double EMin;
gsl_rng * r;

std::ofstream logfile;
std::ofstream integrandfile;
std::ofstream ptfile;
std::ofstream sampfile;
std::ofstream contribfile;
std::ofstream outfile;
std::ofstream outfullfile;
std::ofstream latexfile;
std::ofstream testfile;
int NSamples;
int NTemps;
bool ErrorFlag = false;
std::vector<double> Temp;

bool bEnergyCorrelations=false, bPartialWidthCorrelations=false;

// counters
int PenZeroCount=0, IntegratedCount=0, SubSampledPosCount=0, SampledNegCount=0,
  NANCount=0, InfCount=0, BelowIntLimit=0, IntfNANCount=0, LogZeroCount=0,
  WeirdMuSigmaCount=0;

// Is a string numeric?
bool isNumber(const std::string& str){
	return str.find_first_not_of("0123456789.") == std::string::npos;
}

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

  // Get a line
  std::getline(infile, line);
  std::stringstream ss(line);
  std::string entry;
  // Use c++ stringstream to grab separated elements from the line and
  // put them in a vector of strings
  while( ss >> entry ){
    //    std::cout << entry << "  ";
    entries.push_back(entry);
  }

  // Return to saved position in file
  infile.seekg(place);

  // Return the length of the vector of strings
  return entries.size();
  
}

void readNonResonant(std::ifstream &infile, Reaction &R, int part){

  // First save position
  int place=infile.tellg();

  // Read a bit of data to make sure it's not commented-out
  std::string data;
  infile >> data;
  std::size_t found = data.find("!");
  if (found!=std::string::npos){
    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    R.setNonResonant(0,0,0,0,0,part);
    return;
  }
  // Return to saved position in file
  infile.seekg(place);
  
  int ne = countEntries(infile);
  //std::cout << "There are " << ne << " entries\n";
  if(ne != 5){
    std::cout << "  ERROR: There should be 5 numbers in the non-resonant input lines\n";
    exit(EXIT_FAILURE);
  }

  double s, sp, spp, ds, cutoffe;
  infile >> s >> sp >> spp >> ds >> cutoffe;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R.setNonResonant(s,sp,spp,ds,cutoffe,part);

  if(ds < 0)logfile << "Non-resonant part " << part+1 << " has a factor uncertainty of " 
										<< -ds << std::endl;
  
}

/* 
   This is where Victor need to work on reading all of the "standard" resonances
   1) Start by just reading 5 of them
   2) Then add the ability to add up until the end of the block (check for ****?)
   3) Include error checking. Do the numbers make sense?
*/

void readResonanceBlock(std::ifstream &infile, Reaction &R, bool isUpperLimit){

  double E_cm, dE_cm, wg=0.0, dwg=0.0, Jr,
    G1, dG1, PT1=0.0, DPT1=0.0, G2, dG2, PT2=0.0, DPT2=0.0,
    G3, dG3, PT3=0.0, DPT3=0.0, Exf;
  int i,L1, L2, L3, isBroad;
	bool isECorrelated;

  i=0;
  // Read the number of entries on the first resonance line. This
  // gives us another check on whether it's an upper limit or normal
  // resonance
  int nEnt = countEntries(infile);
  //  std::cout << "Reading resonances:\n";
  //std::cout << "There are " << countEntries(infile) << " entries in this section\n";
  if(!isUpperLimit){
    isUpperLimit = false;
    //std::cout << "Normal resonances\n";
    if(nEnt != 16){
      std::cout << "ERROR! The number of columns in the resonance section is wrong\n";
      std::cout << "       Expect N=16; Got N=" << nEnt << "\n";
    }
  } else {
    isUpperLimit = true;
    //std::cout << "Upper limit resonances\n";
    if(nEnt == 17){
      std::cout << "WARNING: It looks like you're using an old version of RatesMC input without DPT\n";
    } else if (nEnt == 20) {
    } else {
      std::cout << "ERROR! The number of columns in the upper limit resonance section is wrong\n";
      std::cout << "       Expect N=17 or 20; Got N=" << nEnt << "\n";
    }
  }
     
  while(true){

    // First try to read resonance energy to see if it's a real resonance input
    std::string data;
    infile >> data;

    // Look for '***', which signifies the end of resonance input
    std::size_t found;
    found = data.find("***");
    if (found!=std::string::npos){
      logfile << "Found end of ";
      if(isUpperLimit)
				logfile << "upper limit ";
      logfile << "resonances\n" << data << std::endl;
      break;
    }

    // Now look for an exclamation point, which indicates a commented-out resonance
    //logfile << data << "\n";
    found = data.find("!");
    if (found!=std::string::npos){
      if(isUpperLimit)
				logfile << "Upper limit ";
      logfile << "Resonance is commented-out (" << data << ")" << std::endl;
      infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
      continue;
    }
    
    // If this looks like a resonance, read a single line
    E_cm = std::stod(data);
    //    logfile << E_cm << std::endl;;
    if(!isUpperLimit){
			// First check whether the energy uncertainty has a correlation tag
			infile >> data;
			isECorrelated = (data.find("c")!=std::string::npos);
			dE_cm = std::stod(data);
			infile >> wg >> dwg >> Jr 
						 >> G1 >> dG1 >> L1 >> G2 >> dG2 >> L2 >> G3 >> dG3 >> L3
						 >> Exf >> isBroad;
    } else {
      // Upper limit resonances.
      // Old style with no DPT
      if(nEnt == 17){
				// First check whether the energy uncertainty has a correlation tag
				infile >> data;
				isECorrelated = (data.find("c")!=std::string::npos);
				dE_cm = std::stod(data);
				infile >> Jr 
							 >> G1 >> dG1 >> L1 >> PT1 >> G2 >> dG2 >> L2 >> PT2
							 >> G3 >> dG3 >> L3 >> PT3
							 >> Exf >> isBroad;
	
      }
      // New style with DPT
      else if(nEnt == 20){
				// First check whether the energy uncertainty has a correlation tag
				infile >> data;
				isECorrelated = (data.find("c")!=std::string::npos);
				dE_cm = std::stod(data);
				infile >> Jr 
							 >> G1 >> dG1 >> L1 >> PT1 >> DPT1 >> G2 >> dG2 >> L2 >> PT2 >> DPT2
							 >> G3 >> dG3 >> L3 >> PT3 >> DPT3
							 >> Exf >> isBroad;
      }
    }

    // Convert to correct units
    E_cm  *= 1.0e-3;   // keV to MeV
    dE_cm *= 1.0e-3;   // keV to MeV
    wg    *= 1.0e-6;   // eV to MeV
    dwg   *= 1.0e-6;   // eV to MeV
    G1    *= 1.0e-6;   // eV to MeV
    dG1   *= 1.0e-6;   // eV to MeV
    G2    *= 1.0e-6;   // eV to MeV
    dG2   *= 1.0e-6;   // eV to MeV
    G3    *= 1.0e-6;   // eV to MeV
    dG3   *= 1.0e-6;   // eV to MeV

    // if E_cm is negative, G1 is actually unitless, so undo the unit operation above
    if(E_cm < 0){
      G1 *= 1.0e6;
      dG1 *= 1.0e6;
    }
		// If factor uncertainties are input, don't scale 
		if(dG1 < 0.0) dG1 *= 1.0e6;
		if(dG2 < 0.0) dG2 *= 1.0e6;
		if(dG3 < 0.0) dG3 *= 1.0e6;
		if(dwg < 0.0) dwg *= 1.0e6;
		

    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

		// If we care about energy correlations
		if(bEnergyCorrelations && isECorrelated){
			//std::cout << "Energy correlations are enabled!" << std::endl;
			//std::cout << dE_cm << " " << R.smallestdE << " " <<  isZero(R.smallestdE) << "\n";
			if(dE_cm < R.smallestdE || isZero(R.smallestdE))R.smallestdE = dE_cm;
			//			std::cout << R.smallestdE << "\n";
		}
		if(bPartialWidthCorrelations && !isUpperLimit){
			if(!isZero(wg)){
				if(dwg > 0.0){
					if(dwg/wg < R.smallestdwg || isZero(R.smallestdwg))
						R.smallestdwg = dwg/wg;
				} else {
					if( ((-dwg)-1.0) < R.smallestdwg || isZero(R.smallestdwg))
						R.smallestdwg = (-dwg)-1.0;
				}
			}
			
			if(!isZero(G1)){
				if(dG1 > 0.0){
					if(dG1/G1 < R.smallestdG[0] || isZero(R.smallestdG[0]))
						R.smallestdG[0] = dG1/G1;
				} else {
					if( ((-dG1)-1.0) < R.smallestdG[0] || isZero(R.smallestdG[0]))
						R.smallestdG[0] = (-dG1)-1.0;
				}
			}

			if(!isZero(G2)){
				if(dG2 > 0.0){
					if(dG2/G2 < R.smallestdG[1] || isZero(R.smallestdG[1]))
						R.smallestdG[1] = dG2/G2;
				} else {
					if( ((-dG2)-1.0) < R.smallestdG[1] || isZero(R.smallestdG[1]))
						R.smallestdG[1] = (-dG2)-1.0;
				}
			}

			if(!isZero(G3)){
				if(dG3 > 0.0){
					if(dG3/G3 < R.smallestdG[2] || isZero(R.smallestdG[2]))
						R.smallestdG[2] = dG3/G3;
				} else {
					if( ((-dG3)-1.0) < R.smallestdG[2] || isZero(R.smallestdG[2]))
						R.smallestdG[2] = (-dG3)-1.0;
				}
			}
		}
		
		
    // Add this resonance to the list of resonances stored in the
    // reaction
    R.addResonance(i++, E_cm, dE_cm, wg, dwg, Jr,
									 G1, dG1, L1, PT1, DPT1,
									 G2, dG2, L2, PT2, DPT2,
									 G3, dG3, L3, PT3, DPT3,
									 Exf, isBroad, isUpperLimit,isECorrelated);
    
  }
  
}

int ReadInputFile(std::string inputfilename, Reaction *R){

	// Load in the AME Reader
	std::string amefilename = "mass_1.mas20";
	std::string nubasefilename = "nubase_3.mas20";
	AMEReader *ame = new AMEReader(amefilename,nubasefilename);
	
  std::cout << "The input file name is: " << inputfilename << std::endl;

  //int ne;  // for counting entries
  
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
  //R -> getName();

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
  
  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	// Use AME if the input isn't a number
  if(isNumber(dummy)){
		z0 = std::stoi(dummy);
	} else {
		z0 = ame -> readCharge(dummy);
	}

  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	// Use AME if the input isn't a number
  if(isNumber(dummy)){
		z1 = std::stoi(dummy);
	} else {
		z1 = ame -> readCharge(dummy);
	}

	infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	// Use AME if the input isn't a number
  if(isNumber(dummy)){
		z2 = std::stoi(dummy);
	} else {
		z2 = ame -> readCharge(dummy);
	}

	R -> setCharges(z0,z1,z2);

	// M0 (Projectile mass)
  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	// Use AME if the input isn't a number
  if(isNumber(dummy)){
		m0 = std::stod(dummy);
	} else {
		m0 = ame -> readMass(dummy);
		m0 = atomicToNuclear(m0, z0);
		std::cout << "= " << m0 << " (nuclear)\n";
	}
	// M1 (Target mass)
	infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		m1 = std::stod(dummy);
	} else {
		m1 = ame -> readMass(dummy);
		m1 = atomicToNuclear(m1, z1);
		std::cout << "= " << m1 << " (nuclear)\n";
	}
	// M2 (Ejectile mass)
	infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		m2 = std::stod(dummy);
	} else {
		m2 = ame -> readMass(dummy);
		m2 = atomicToNuclear(m2, z2);
		std::cout << "= " << m2 << " (nuclear)\n";
	}
  R -> setMasses(m0,m1,m2);

  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		j0 = std::stod(dummy);
	} else {
		j0 = ame -> readSpin(dummy);
	}

	infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		j1 = std::stod(dummy);
	} else {
		j1 = ame -> readSpin(dummy);
	}

  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		j2 = std::stod(dummy);
	} else {
		j2 = ame -> readSpin(dummy);
	}

	R -> setSpins(j0,j1,j2);
	//  logfile << "   " << j0 << "    " << j1 << "    " << j2 << "\n";

  // Entrance and exit particle separation energies
	infile >> dummy;
	infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		Qin = std::stod(dummy);
	} else if (dummy == "ame" || dummy == "AME"){
		int acompound = (int)round(m0+m1);
		int zcompound = z0+z1;
		double mcompound = ame -> readMassFromAandZ(acompound, zcompound);
		mcompound = atomicToNuclear(mcompound, zcompound);
		Qin = ((m0 + m1) - mcompound)*AMU*1000.0;   // To get into keV
	} else {
		std::cout << "ERROR: Enter a number, 'ame', or 'AME'"
							<< " for the entrance particle separation energy\n";
		exit(EXIT_FAILURE);
	}
	//	std::cout << "Qin = " << Qin << "\n";
	
	infile >> dummy;
	infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(isNumber(dummy)){
		Qout = std::stod(dummy);
	} else if (dummy == "ame" || dummy == "AME"){
		// First find the compund nucleus
		int acompound = (int)round(m0+m1);
		int zcompound = z0+z1;
		double mcompound = ame -> readMassFromAandZ(acompound, zcompound);
		mcompound = atomicToNuclear(mcompound, zcompound);

		// Then find the residual nucleus
		int aresidual = (int)round(m0+m1-m2);
		int zresidual = z0+z1-z2;
		double mresidual = ame -> readMassFromAandZ(aresidual, zresidual);
		mresidual = atomicToNuclear(mresidual, zresidual);

		Qout = ((m2 + mresidual) - mcompound)*AMU*1000.0;   // To get into keV
	} else {
		std::cout << "ERROR: Enter a number, 'ame', or 'AME'"
							<< " for the exit particle separation energy\n";
		exit(EXIT_FAILURE);
	}
	//	std::cout << "Qout = " << Qout << "\n";

	
	//  Qout = readDouble(infile);
  R -> setSeparationEnergies(Qin/1000.0, Qout/1000.0);

  // The R0 radius
  R0 = readDouble(infile);
  R -> setR0(R0);

  // Index of the gamma-ray channel
  gindex = readInt(infile);
  R -> setGammaIndex(gindex-1);

  // Ignore a line
  skipLines(infile, 1);

  // Computation control block
  EMin = readDouble(infile)/1000.;
  
  NSamples = readInt(infile);
  NTemps = readInt(infile);

  // Only read correlations input if it exists
  int place=infile.tellg();
  infile >> dummy;
  std::size_t found = dummy.find("***");
  if(found!=std::string::npos){
    std::cout << "WARNING: It looks like you are using an old RatesMC.in.\n" <<
      "         I'll assume you don't care about correlations!\n" << std::endl;
    bPartialWidthCorrelations = false;
    bEnergyCorrelations = false;

  } else {
    infile.seekg(place);
    int itmp;
    infile >> itmp;
    bPartialWidthCorrelations = (bool)itmp;
    infile >> itmp;
    bEnergyCorrelations = (bool)itmp;
    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    skipLines(infile, 1);
  }
  
  // Skip 2 lines
  skipLines(infile, 2);

  // Non-resonant line 1
  readNonResonant(infile, *R, 0);
  
  // Non-resonant line 2
  readNonResonant(infile, *R, 1);

  // Skip 5 lines
  skipLines(infile, 5);

  // Read all of the resonances 
  readResonanceBlock(infile, *R, false);

  // Skip 4 lines
  skipLines(infile, 4);
  // Read the upper limit resonances
  readResonanceBlock(infile, *R, true);

  logfile << std::endl;


	delete(ame);
	
  return 0;
}

//----------------------------------------------------------------------
// Define the temperature array
void defineTemperatures(){

  
   std::vector<double> defaultT{0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
	 														 0.01,0.011,0.012,0.013,0.014,0.015,
	 														 0.016,0.018,0.020,0.025,0.03,0.04,
	 														 0.05,0.06,0.07,0.08,0.09,0.1,0.11,
	 														 0.12,0.13,0.14,0.15,0.16,0.18,0.20,
	 														 0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,
	 														 0.8,0.9,1.0,1.25,1.5,1.75,2,2.5,3,
	 														 3.5,4,5,6,7,8,9,10};
  
  
  //std::vector<double> defaultT{0.001,0.01,0.05,0.1};
  Temp = defaultT;

  logfile << "--------------------------------------------------\n";
  logfile << "There are " << Temp.size() << " temperatures:\n";
  for(std::size_t iT=0; iT<Temp.size(); iT++){
    logfile << std::setw(5) << Temp[iT] << "  ";
    if((iT+1)%7 == 0 && iT>0)logfile << "\n";
  }
  logfile << "\n";
  logfile << "--------------------------------------------------";
  logfile << std::endl;

}

//----------------------------------------------------------------------
void logNormalize(double mean, double sd, double& mu, double& sigma){

  // Varience is standard deviation squared
  double var = gsl_pow_2(sd);

  if(mean == 0.0){
    mu = 0.0;
    sigma = 0.0;
    return;
  }

  mu = gsl_sf_log(mean) - 0.5*gsl_sf_log(1+(var/gsl_pow_2(mean)));
  sigma = sqrt(gsl_sf_log(1+(var/gsl_pow_2(mean))));

  // If mean and var are very small, some computers will return nan for these ratios
  // so check for it
  if(gsl_isnan(mu) || gsl_isnan(sigma)){
    mu = 0.0;
    sigma = 0.0;
  }

  return;


}


//----------------------------------------------------------------------
// Penetration factor calculation
double PenFactor(double E, double L, double Mass0, double Mass1,
								 int Charge0, int Charge1, double R){

  //  std::cout << E << " " << L << " " << Mass0 << " " << Mass1 << " " << Charge0
  //	    << " " << Charge1 << " " << R << "\n";
  
  // Turn off the GSL error handler, which aborts the program
  // if G or F go out of range.
  gsl_set_error_handler_off();

  gsl_sf_result F,Fp,G,Gp;
  double exp_F,exp_G,P;
  double mue = Mass0*Mass1/(Mass0+Mass1);

  //  double R = R0*(pow(Mass0,(1./3.))+pow(Mass1,(1./3.)));
  double rho = 0.218735097*R*sqrt(mue*E);

  double eta = 0.15748927*(double)Charge0*(double)Charge1*sqrt(mue/E);

  int status = gsl_sf_coulomb_wave_FG_e (eta, rho, L, 0, &F, &Fp, &G,
																				 &Gp, &exp_F, &exp_G);

  // Check to see if this failed. If it's out of range, set P=0,
  // if not, print an error and exit
  if(status){
    if(status == GSL_EOVRFLW){
      ErrorFlag = true;
      PenZeroCount++;
      return 0.0;
    } else {
      std::cout << "\nERROR: Something went wrong in coulomb wavefunction!" <<
				"\n\tGSL Error: " << gsl_strerror (status)<< std::endl;
      std::cout << "The Energy was " << E*1e3 << " keV." << std::endl;
      abort();
    }
  }
  gsl_set_error_handler (NULL);


  // Just in case there is an overflow, multiply by exponential
  // (See GSL documentation for more info)
  double F_l = F.val*exp(exp_F);
  double G_l = G.val*exp(exp_G);

  P = rho/(gsl_pow_2(F_l) + gsl_pow_2(G_l));

  return P;
}

//----------------------------------------------------------------------
// Transpose a 2D vector of doubles
void transpose(std::vector<std::vector<double> > &b){

  if(b.size() == 0)
    return;

  std::vector<std::vector<double> > trans_vec(b[0].size(), std::vector<double>());

  for (std::size_t i = 0; i < b.size(); i++){
    for (size_t j = 0; j < b[i].size(); j++) {
      trans_vec[j].push_back(b[i][j]);
    }
  }

  b = trans_vec;    // <--- reassign here
  
}

//----------------------------------------------------------------------
void writeOutputFileHeaders(Reaction *R){

  outfile << R->getName() << std::endl;
  outfile << "Calculated with RatesMC " << VersionNumber << " (" <<
    VersionDate << ")" << std::endl;
  outfile << "Samples = " << NSamples << std::endl;
  outfile << " T9      RRate_low       Median Rate" << 
    "     RRate_high     f.u."  << std::endl;
  
  outfullfile << R->getName() << std::endl;
  outfullfile << "Calculated with RatesMC " << VersionNumber << " (" <<
    VersionDate << ")" << std::endl;
  outfullfile << "Samples = " << NSamples << std::endl;
  outfullfile << " T9      RRate_2low      RRate_low       Classical Rate  Median Rate" << 
    "     Mean Rate       RRate_high     RRate_2high     Log-Normal mu" <<
    "      Log-Normal sigma A-D Statistic" << std::endl;

  sampfile << R->getName() << std::endl;
  sampfile << "Calculated with RatesMC " << VersionNumber << " (" <<
    VersionDate << ")" << std::endl;
  
}

//----------------------------------------------------------------------
void writeContributions(std::vector<std::vector<double> > Contributions, double Temperature){

  //std::cout << "At the end we have\n";
	//std::cout << "RateSample of length: " << RateSample.size() << "\n";
	transpose(Contributions);
	//std::cout << "Contributions of length: " << Contributions.size() << " X "
	//	      << Contributions[0].size() << "\n";
	//std::cout << std::endl;
    
	// From the contributions, calculate Low, median, and high
	// contribution for each resonance
	// Sort, and find quantiles for each rate contribution
	std::vector<double> LowCont;
	std::vector<double> HighCont;
	std::vector<double> MedianCont;
	for(std::size_t j=0;j<Contributions.size();j++){

		// Remove all nan values
		Contributions[j].erase(std::remove_if(std::begin(Contributions[j]),
																					std::end(Contributions[j]),
																					[](const auto& value) { return isnan(value); }),
													 std::end(Contributions[j]));
			
		// Sort the Contributions for each resonance
		std::sort (Contributions[j].begin(), Contributions[j].end());
			
		// Now convert this into a gsl vector
		gsl_vector_const_view gsl_v =
			gsl_vector_const_view_array( &Contributions[j][0], Contributions[j].size() );
      
		LowCont.push_back(gsl_stats_quantile_from_sorted_data(gsl_v.vector.data,1, Contributions[j].size(),0.16));
		HighCont.push_back(gsl_stats_quantile_from_sorted_data(gsl_v.vector.data,1, Contributions[j].size(),0.84));
		MedianCont.push_back(gsl_stats_median_from_sorted_data(gsl_v.vector.data,1,Contributions[j].size()));
	}

    
	// for(int i=0; i<LowCont.size(); i++){
	//   std::cout << "i = " << i << "  LowCont = " << LowCont[i] <<  "  MedianCont = " << MedianCont[i]
	// << "  HighCont = " << HighCont[i] << "\n";
	// }
    
    
	// Write to file
	contribfile.precision(3);
	contribfile.width(5);
	contribfile << Temperature << std::scientific << "   " ;
	for(std::size_t i=0; i<MedianCont.size(); i++){
		if(LowCont[i] < 1.0e-99)LowCont[i]=0.0;
		if(MedianCont[i] < 1.0e-99)MedianCont[i]=0.0;
		if(HighCont[i] < 1.0e-99)HighCont[i]=0.0;
		//if(j<2 || UseInRate[j-2] || j>=NRes+2){
		contribfile << LowCont[i] << "  " << MedianCont[i] << "  " <<
			HighCont[i] << "  ";
	}
    
	contribfile << std::endl;
	contribfile.unsetf(std::ios_base::scientific);
    

  
}

//----------------------------------------------------------------------
// Write the rate to output and LaTeX file
void writeRates(std::vector<double> Rates, double ARate, double Temperature){

  // If there are enough samples, calculate the mean, variance,
  // log-mean and log-variance and bin the rates
  // Turn off the gsl error handler in case we have a negative number here

  double MeanRate, RateMu, RateSigma, fu;
  // The 1- and 2-sigma high and low rates
  double Low2Rate, LowRate, MedianRate, HighRate, High2Rate;
  // Parentheses characters
  char Parenth[2] = {' ',' '};

  
  gsl_set_error_handler_off();
	std::vector<double> logRates;
  gsl_sf_result LogResult;
  if(NSamples > 2){
    for(int s=0;s<NSamples;s++){
      //summed[k] *= RateFactor;
      int status = gsl_sf_log_e(Rates[s],&LogResult);
      // Check to make sure log didn't error
      if(status){
				ErrorFlag = true;
				LogZeroCount++;
				logRates.push_back(std::numeric_limits<double>::quiet_NaN());
      } else {
				logRates.push_back(LogResult.val);
				//				std::cout << lograte[s] << "\n";
      }
    }
    /*    
		// No longer need to bin rates
		//      BinRates(summed,Temp[i],i,distfileName);
		for(int k=0;k<NSamples;k++){
		sampfile << summed[k] << endl;
		summed[k] /= RateFactor;
		}
		sampfile << endl;
    */
  }

	// Remove all nan rates
	Rates.erase(std::remove_if(std::begin(Rates),
														 std::end(Rates),
														 [](const auto& value) { return isnan(value); }),
							std::end(Rates));
	logRates.erase(std::remove_if(std::begin(logRates),
														 std::end(logRates),
														 [](const auto& value) { return isnan(value); }),
							std::end(logRates));
	
	
  // Now, find the uncertainties. Before finding quantiles, need to sort
	std::sort(Rates.begin(), Rates.end());
  // Now convert the rates into a gsl vector
  gsl_vector_const_view gsl_Rates =
    gsl_vector_const_view_array( &Rates[0], Rates.size() );
  gsl_vector_const_view gsl_logRates =
    gsl_vector_const_view_array( &logRates[0], logRates.size() );

	
  // Means and variances
  MeanRate = gsl_stats_mean(gsl_Rates.vector.data,1,Rates.size());
  RateMu = gsl_stats_mean(gsl_logRates.vector.data,1,logRates.size());
  RateSigma = sqrt(gsl_stats_variance(gsl_logRates.vector.data,1,logRates.size()));
  fu = gsl_sf_exp(RateSigma);
    
  // Uncertainties calculated from quantiles
  Low2Rate =  gsl_stats_quantile_from_sorted_data (gsl_Rates.vector.data,1,Rates.size(),0.025);
  LowRate = gsl_stats_quantile_from_sorted_data   (gsl_Rates.vector.data,1,Rates.size(),0.16);
  MedianRate =  gsl_stats_median_from_sorted_data (gsl_Rates.vector.data,1,Rates.size());
  HighRate =  gsl_stats_quantile_from_sorted_data (gsl_Rates.vector.data,1,Rates.size(),0.84); 
  High2Rate =  gsl_stats_quantile_from_sorted_data(gsl_Rates.vector.data,1,Rates.size(),0.975); 
  
  // Turn off the gsl error handler in case we have a negative number here
  gsl_set_error_handler_off();
  gsl_sf_result LogMean;

  int status = gsl_sf_log_e(MeanRate,&LogMean);

  // if log calculation failed, set mu and sigma to zero and print error
  if(status){
    ErrorFlag = true;
    std::cout << "\nERROR: Something went wrong in logarithm, negative or zero rate?!\n"
							<< std::endl;
    RateMu = 0.0;
    RateSigma = 0.0;
  } else {
    // Calculate mu and sigma 
    // Check to see if mu is close to ln(median)
    double mu_err = (RateMu-gsl_sf_log(MedianRate))/
      gsl_sf_log(MedianRate);

    // sigma error using high and medium rates
    double sigma_err =
      (RateSigma - gsl_sf_log(HighRate/MedianRate))/
      gsl_sf_log(HighRate/MedianRate);
    
    
    if(fabs(mu_err) > 0.03 || fabs(sigma_err)>0.2){
      ErrorFlag = true;
      WeirdMuSigmaCount++;
      
      Parenth[0] = ' ';
      Parenth[1] = ' ';
    }

  }

  
  if(LogZeroCount > 0){
    logfile << "\tWARNING: The Rate was zero, creating \n\t\tproblems in ln " <<
      LogZeroCount << " times" << std::endl;
  }

  // Check validity of lognormal. Only for more than 1 sample
  double AndDar_Asqrd=0.0;
  if(NSamples > 1){
    AndDar_Asqrd = CalcAD(Rates,RateMu,RateSigma);
  } else {
    AndDar_Asqrd = 0.0;
  }

  //  std::cout << "AD = " << AndDar_Asqrd << std::endl;
  // Print stuff
  std::cout << "\n";
  std::cout << "For T9 = " << Temperature << " the median rate = ";
  std::cout.precision(3);
  std::cout << std::scientific << MedianRate << std::endl;
  std::cout.unsetf(std::ios_base::scientific);


  if(LowRate < 1.0e-99)LowRate=0.0;
  if(ARate < 1.0e-99)ARate=0.0;
  if(Low2Rate < 1.0e-99)Low2Rate=0.0;
  if(High2Rate < 1.0e-99)High2Rate=0.0;
  if(MedianRate < 1.0e-99)MedianRate=0.0;
  if(MeanRate < 1.0e-99)MeanRate=0.0;
  if(HighRate < 1.0e-99)HighRate=0.0;
    
  char buffer[200];
  // RatesMC.out output
  sprintf(buffer,"%6.3f  %10.3e      %10.3e      %10.3e       %10.3e",
					Temperature,LowRate,MedianRate,HighRate,fu);
  outfile << buffer << std::endl;
  //  outfile << std::endl;
  
  // RatesMC.full output
  sprintf(buffer,"%6.3f  %10.3e      %10.3e      %10.3e      %10.3e      %10.3e      %10.3e     %10.3e      %c% 10.3e%c      %c% 10.3e%c       %9.3e",
					Temperature,Low2Rate, LowRate, ARate, MedianRate, MeanRate, HighRate,
					High2Rate, Parenth[0],RateMu,Parenth[1],Parenth[0],RateSigma,Parenth[1],AndDar_Asqrd);
  outfullfile << buffer << std::endl;
  //  outfullfile << std::endl;

  WriteLatex2(Temperature, LowRate, MedianRate, HighRate, RateSigma);

}

//----------------------------------------------------------------------
// Function to write a latex table
void WriteLatex2(double Temperature, double LowRate, double MedianRate, double HighRate,
								 double RateSigma){

  double low_x,low_f,median_x,median_f,high_x,high_f,fu;

  // Write a test rate to file
  //   ofile << reactionname << endl;
  //   ofile << "Samples = " << NSamples << endl;
  //   ofile << "T9 & Min Rate & Mean Rate & Max Rate & mu & sigma \\\\" << endl;

  // set the precision
  latexfile.precision(3);
  latexfile.setf(std::ios::fixed,std::ios::floatfield);
  latexfile << std::setfill('0');
  latexfile << std::setiosflags(std::ios::internal);
  //  for(int i=0;i<NTemps;i++){

  // Calculate all of the exponent stuff
  low_f = floor(log(LowRate)/log(10.0));
  low_x = LowRate*pow(10.0,-low_f);
  median_f = floor(log(MedianRate)/log(10.0));
  median_x = MedianRate*pow(10.0,-median_f);
  high_f = floor(log(HighRate)/log(10.0));
  high_x = HighRate*pow(10.0,-high_f);
  //mu_sign = (RateMu>0)-(RateMu<0);
  //mu_f = floor(log(fabs(RateMu))/log(10.0));
  //mu_x = RateMu*pow(10.0,-mu_f);
  //sigma_f = floor(log(RateSigma)/log(10.0));
  //sigma_x = RateSigma*pow(10.0,-sigma_f);
  //AD_f = floor(log(AD)/log(10.0));
  //AD_x = AD*pow(10.0,-AD_f);
  fu = exp(RateSigma);

  //int prec = 3 - floor(log10(fu));
    
  latexfile << Temperature << " & ";
	if(isnan(median_x)){
		low_x = 0.0;
		low_f = 0.0;
		median_x = 0.0;
		median_f = 0.0;
		high_x = 0.0;
		high_f = 0.0;
		fu = 1.0;
	}
		
	latexfile << std::setprecision(3) << low_x <<"E"
						<< std::setprecision(0) << std::setw(3)
						<< std::setiosflags(std::ios::showpos)
						<< low_f << " & " << std::setprecision(3)
						<< std::resetiosflags(std::ios::showpos)
						<< median_x << "E" << std::setw(3)
						<< std::setiosflags(std::ios::showpos)
						<< std::setprecision(0) << median_f << "& \n" << std::setprecision(3)
						<< std::resetiosflags(std::ios::showpos)
						<< "      " << high_x << "E" << std::setw(3)
						<< std::setiosflags(std::ios::showpos)
						<< std::setprecision(0) << high_f << " & " << std::setprecision(3)
						<< std::resetiosflags(std::ios::showpos)
						<< fu << " \\\\ " << std::setprecision(3)
						<< std::resetiosflags(std::ios::showpos) << std::endl;
	
  
  //  latexfile << "\n";
  //latexfile.close();

  return;
}

//----------------------------------------------------------------------
void writeRateSamples(std::vector<double> RateSample, double Temp){

  sampfile << "T9 = " << Temp << std::endl;
  sampfile << "Samples = " << RateSample.size() << std::endl;

  for(double rate : RateSample)
    sampfile << rate << "\n";

  sampfile << std::endl;
  return;
}

//----------------------------------------------------------------------
void summarizeErrors(double Temp){

  // For this temperature, let the user know if nothing went wrong!
  if(!ErrorFlag) {
    logfile << "\tNo Errors Occured!" << std::endl;
  } else {

    if(SampledNegCount>0){
      logfile << " A positive resonance had a negative energy "
							<< SampledNegCount << " times" << std::endl;
      SampledNegCount = 0;
    }
    if(SubSampledPosCount>0){
      logfile << " A negative resonance had a positive energy "
							<< SubSampledPosCount << " times" << std::endl;
      SubSampledPosCount = 0;
    }
    if(IntegratedCount>0){
      logfile << " A positive narrow resonance had a negative energy "
							<< IntegratedCount << " times" << std::endl;
      IntegratedCount = 0;
    }
    if(NANCount>0){
      logfile << " An integration error occured "
							<< NANCount << " times" << std::endl;
      NANCount = 0;
    }
    if(InfCount>0){
      logfile << " An integration reached infinity "
							<< InfCount << " times" << std::endl;
      InfCount = 0;
    }
    ErrorFlag = false;
  }

  logfile << "----------------------------------------\n";
  
}

//----------------------------------------------------------------------
#include <gsl/gsl_sf_erf.h>
//---------------------------------------------------------------
// Check how well the lognormal fits
double CalcAD(std::vector<double> Rates,double Mu,double Sigma){

  // For now, do a Kolmogorov-Smirnov Goodness-of-Fit Test
  double KS = 0.0;
  double CumLognorm[Rates.size()+1];
  double X,maxarray[3];
  double AndDar_S=0.0, AndDar_Asqrd=0.0;
  int NanADCount = 0;

  // For the anderson-Darling test, convert mu and sigma to normal parameters
  // and take the log of each point
  //double DistMean = gsl_sf_log(gsl_sf_exp(Mu + 0.5*gsl_pow_2(Sigma)));
  //double DistVar = gsl_sf_log((gsl_sf_exp(gsl_pow_2(Sigma))-1.0)*
  //		      gsl_sf_exp(2.0*Mu + gsl_pow_2(Sigma)));

  CumLognorm[0]=0.0;
  for(std::size_t i=1;i<(Rates.size()+1);i++){
    X = Rates[i-1];
    // if X=0, can't take log, so check
    if(X==0.0){
      CumLognorm[i] = CumLognorm[i-1];
    } else {
      CumLognorm[i] = CumLognorm[i-1]+((1.0/(X*Sigma*sqrt(2.0*M_PI)))*
																			 gsl_sf_exp(-gsl_pow_2(gsl_sf_log(X)-Mu)/
																									(2.0*gsl_pow_2(Sigma))));
    }
  }
  double norm = 1/CumLognorm[Rates.size()];
  double AndDar_S_tmp;

  for(std::size_t i=1;i<(Rates.size()+1);i++){
    X = Rates[i-1];

    // check to make sure X isn't zero
    if(X==0.0){
      AndDar_S_tmp = 0.0;
    } else {
      // also calculate the Anderson-Darling test
      double F = 0.5*(gsl_sf_erf((gsl_sf_log(Rates[i-1])-Mu)/
																 (sqrt(2.0)*Sigma)) + 1.0);
      double OneMinusF = 1.0-
				0.5*(gsl_sf_erf((gsl_sf_log(Rates[Rates.size()-i])-Mu)/
												(sqrt(2.0)*Sigma)) + 1.0);

      AndDar_S_tmp =
				((2.0*i-1.0)/Rates.size())*(gsl_sf_log(F)+gsl_sf_log(OneMinusF));
    }

    // check for NaN (if F is too small)
    if(!isnan(AndDar_S_tmp)){
      AndDar_S += AndDar_S_tmp;
    } else {
      NanADCount++;
    }

    maxarray[0]=KS;
    maxarray[1]=norm*CumLognorm[i] - double(i-1.0)/Rates.size();
    maxarray[2]=double(i)/Rates.size() - norm*CumLognorm[i];
    KS = gsl_stats_max(maxarray,1,3);
  }

  AndDar_Asqrd = (-(double)Rates.size()-AndDar_S)*(1+4.0/(double)Rates.size()+25.0/(double)gsl_pow_2(Rates.size()));

  if(NanADCount > 0){
    logfile << "\tWARNING: Anderson-Darling test cannot \n\t\tbe calculated for " << NanADCount <<
      " sample(s).\n" << std::endl;
    ErrorFlag=true;
  }

  return AndDar_Asqrd;
}

//----------------------------------------------------------------------
// Check for zero
bool isZero(double x){
  const double epsilon = 1e-5;
  return fabs(x-0.) <= epsilon*fabs(x);
}

//----------------------------------------------------------------------
// Convert atomic to nuclear mass
double atomicToNuclear(double A, double Z) {
	return A - Z*ElectronMass;
}

//----------------------------------------------------------------------
// Check if a file exists
inline bool fileExists (const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

//----------------------------------------------------------------------
// Setup the random sampler
void setupRandom(){
  
	// set-up the GSL random sampler
  const gsl_rng_type * T;
  /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  // Get a seed from either /dev/random, or the clock
  unsigned long int seed = random_seed();
  // Set the seed in the RNG
  gsl_rng_set(r,seed);

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Function to find the random seed
// Courtesy of Robert G. Brown

#include <sys/time.h>

unsigned long int random_seed()
{

  unsigned int seed;
  int iret=0;
  struct timeval tv;
  FILE *devrandom;

  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
    //   if(verbose == D_SEED) printf("Got seed %u from gettimeofday()\n",seed);
  } else {
    iret = fread(&seed,sizeof(seed),1,devrandom);
    //   if(verbose == D_SEED) printf("Got seed %u from /dev/random\n",seed);
    fclose(devrandom);
  }
  if(iret < 0)std::cout << "WARNING! There could be a problem with the random numbers!\n";
  return(seed);

}

//----------------------------------------------------------------------
// Error handler
void RatesMC_Error_Handler(const char * reason,
													 const char * file,
													 int line,
													 int gsl_errno){};
