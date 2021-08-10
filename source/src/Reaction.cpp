/* ======================================================================
   Reaction.cpp
   Author: R. Longland
   Date: 2018-11-13
   
   Description: Contains all of the reaction specific stuff
   ======================================================================
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "stdio.h"

#include <gsl/gsl_randist.h>

#include "Resonance.h"
#include "Reaction.h"
#include "Utilities.h"

using std::cout;
using std::endl;

Reaction::Reaction(){}
Reaction::~Reaction(){}

void Reaction::setNonResonant(double s, double sp, double spp, double ds, double cutoffe, int part){
  S[part] = s;
  Sp[part] = sp;
  Spp[part] = spp;
  dS[part] = ds;
  CutoffE[part] = cutoffe;

  smallestdE = 0.0;
  smallestdwg = 0.0;
}

void Reaction::addResonance(int i, double E_cm, double dE_cm, double wg, double dwg, double Jr,
			    double G1, double dG1, int L1, double PT1, double dPT1,
			    double G2, double dG2, int L2, double PT2, double dPT2,
			    double G3, double dG3, int L3, double PT3, double dPT3,
			    double Exf, bool bInt, bool bUpperLimit){

  double G[3] = {G1, G2, G3};
  double dG[3] = {dG1, dG2, dG3};
  int L[3] = {L1, L2, L3};
  double PT[3] = {PT1, PT2, PT3};
  double dPT[3] = {dPT1, dPT2, dPT3};
  
  // Make a resonance
  Resonance Res(i, E_cm, dE_cm, wg, dwg, Jr,
		G, dG, L, PT, dPT, Exf, bInt, bUpperLimit);
  //Res.print();
  // Add that resonance to the list of resonances
  Resonances.push_back(Res);
}

void Reaction::getName(){
  cout << "The reaction name is: " << Name << "\n";
}

void Reaction::printReaction(){

  cout << "--------------------------------------------------" << "\n";
  cout << "     This is reaction: " << Name << "\n";
  cout << "   Z0        Z1       Z2" << "\n";
  printf( "   %2d        %2d       %2d\n",Z0,Z1,Z2);
  cout << "   M0        M1       M2" << "\n";
  printf( "%5.2f     %5.2f    %5.2f\n",M0,M1,M2);
  cout << "   S_entrance = " << Q << "\n";
  cout << "   S_exit     = " << Qexit << "\n";
  cout << "   The gamma ray channel is channel " << Gamma_index << "\n";
  cout << "--------------------------------------------------" << "\n";

  cout << " Direct Capture part     \n";
  cout << "          S0       S'       S''    dS   CutoffE\n";
  cout << " Part 1: " << S[0] << "  " << Sp[0] << "  " << Spp[0] << "  " << dS[0] << "  " << CutoffE[0] << "\n";
  cout << " Part 2: " << S[1] << "  " << Sp[1] << "  " << Spp[1] << "  " << dS[1] << "  " << CutoffE[1] << "\n";
  cout << "\n";

  cout << "--------------------------------------------------" << "\n";
  cout << " Resonances:\n";
  // Loop through all regular resonances
  std::vector<Resonance>::iterator res;
  for(res = Resonances.begin(); res < Resonances.end(); res++){
    res->print();
  }

}

void Reaction::writeReaction(){

  logfile << "--------------------------------------------------" << "\n";
  logfile << "     This is reaction: " << Name << "\n";
  logfile << "   Z0        Z1       Z2" << "\n";
  logfile << std::setw(5) << Z0 << "     " << std::setw(5) << Z1 << "    " << std::setw(5) << Z2 << "\n";
  logfile << "   M0        M1       M2" << "\n";
  logfile << M0 << "   " << M1 << "   " << M2 << "\n";
  logfile << "   S_entrance = " << Q << "\n";
  logfile << "   S_exit     = " << Qexit << "\n";
  logfile << "   The gamma ray channel is channel " << Gamma_index << "\n";
  logfile << "--------------------------------------------------" << "\n";

  logfile << " Direct Capture part     \n";
  logfile << "          S0       S'       S''    dS   CutoffE\n";
  logfile << " Part 1: " << S[0] << "  " << Sp[0] << "  " << Spp[0] << "  " << dS[0] << "  " << CutoffE[0] << "\n";
  logfile << " Part 2: " << S[1] << "  " << Sp[1] << "  " << Spp[1] << "  " << dS[1] << "  " << CutoffE[1] << "\n";
  logfile << "\n";

  logfile << "--------------------------------------------------" << "\n";
  logfile << " Resonances:\n";
  // Loop through all regular resonances
  std::vector<Resonance>::iterator res;
  for(res = Resonances.begin(); res < Resonances.end(); res++){
    res->write();
  }

}

//----------------------------------------------------------------------
// Calculate the direct capture, non-resonant part of the reaction rate
// NOT WORKING!
double Reaction::calcNonResonant(){

  /*
  double mue = M0*M1/(M0+M1);
  double cutoff_T = 19.92*pow(CutoffE[j],1.5)/sqrt(pow(Z0*Z1,2.)*mue);

  double ADRate,mu,sigma;
  // e designates exact numbers
  double C1e,C2e,C3e,C4e,C5e,C6e,C7e;

  C1e = 7.8324e9*pow((pow(Z0*Z1,2.)*mue),1./6.) *S[j] / sqrt(mue);
  C2e = 4.2475*pow((pow(Z0*Z1,2.)*mue),1./3.);
  C3e = 9.810e-2*pow((pow(Z0*Z1,2.)*mue),-1./3.);
  C4e = 0.1220*(Sp[j]/S[j])*pow((pow(Z0*Z1,2.)*mue),1./3.);
  C5e = 8.377e-2*(Sp[j]/S[j]);
  C6e = 7.442e-3*(Spp[j]/S[j])*pow((pow(Z0*Z1,2.)*mue),2./3.);
  C7e = 1.299e-2*(Spp[j]/S[j])*pow((pow(Z0*Z1,2.)*mue),1./3.);
  */

  /*cout << C1e << "\t" << C2e << "\t" << C3e << "\t" << C4e << "\t" <<
    C5e << "\t" << C6e << "\t" << C7e << endl;
    cout << "Cutoff T = " << cutoff_T << endl;*/

  // New way of doing things: calculate the rate and then sample from
  // the end.

  /*
  // Turn off the error handler in case exp returns zero
  gsl_set_error_handler_off();
  ADRate = (C1e/pow(Temp,2./3.))*gsl_sf_exp(-C2e/pow(Temp,1./3.)-
					    pow(Temp/cutoff_T,2.))*
    (1 + C3e*pow(Temp,1./3.) + C4e*pow(Temp,2./3.) + C5e*Temp +
     C6e*pow(Temp,4./3.) + C7e*pow(Temp,5./3.));

  gsl_set_error_handler(NULL); // Turn error handler back on again
  */
  
  /*cout << (C1e/pow(Temp,2./3.))*exp(-C2e/pow(Temp,1./3.)-
    pow(Temp/cutoff_T,2.)) << endl;
    cout << C3e*pow(Temp,1./3.) << endl;
    cout << C4e*pow(Temp,2./3.) << endl;
    cout <<  C5e*Temp << endl;

    cout << "ADRate = " << ADRate << endl;
    cout << "About to Lognormalise" << endl;*/

  /*
  // If ADRate is negative, set to zero, this is unphysical
  if(ADRate < 0.){
    ErrorFlag = true;
    logfile << "\tWARNING: The non-resonant part caused a negative rate, \n\t\tsetting to zero." << endl;
    ADRate = 0.0;
    for(int i=0;i<NSamples;i++){
      DRate[i]=0.0;
    }
  } else {
    LogNormalise(ADRate,dS[j]*ADRate,mu,sigma);
    // if mu and sigma return as zero, the rate is zero (very small)
    if((mu==0. && sigma==0.)){
      ADRate=0.0;
      for(int i=0;i<NSamples;i++){
	DRate[i]=0.0;
      }
    } else {
      for(int i=0;i<NSamples;i++){
	DRate[i] = gsl_ran_lognormal(r,mu,sigma)/(1.5399e11/pow(mue*Temp,1.5));
	//cout << DRate[i] << endl;
      }
    }
  }


  ADRate = ADRate/(1.5399e11/pow(mue*Temp,1.5));
  */
  return 0.0;
}

//----------------------------------------------------------------------
// Prepare the Monte Carlo samples. 
void Reaction::prepareSamples(){

  cout << "Preparing " << NSamples << " samples\n";

  // The reference samples used for energies, gamma widths, and
  // resonances. There are three sets for each partial width. The
  // first ones are recycled for resonance strengths with the
  // assumption that those are correlated with entrance channel
  // partial widths
  //
  // These are stored by [row][column] where each row is a sample
  std::vector<double> row;
  row.resize(4);
  for(int s=0;s<NSamples;s++){
    for(int j=0; j<4; j++){
      row[j] = gsl_ran_gaussian(r,1.0);
    }
    Ref_sample.push_back(row);
  }
  /*
  for(int s=0; s<NSamples; s++){
    for(int j=0; j<3; j++){
      cout << Ref_sample[s][j] << "  ";
    }
    cout << "\n";
  }
  */
  // For each resonance, go through and calculate all random samples
  for(Resonance Res : Resonances){
    Res.makeSamples(Ref_sample, smallestdE, smallestdwg);
  }

  
}
