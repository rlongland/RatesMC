/* ======================================================================
   Resonance.cpp
   Author: R. Longland
   Date: 2021-01-27
   
   Description: Contains all of the resonance specific stuff
   ======================================================================
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include "stdio.h"

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>

#include "Resonance.h"
#include "Utilities.h"

using std::cout;
using std::endl;

Resonance::Resonance(int index, double E_cm, double dE_cm, double wg, double dwg, double Jr,
		     double G[3], double dG[3], int L[3], double PT[3], double dPT[3],
		     double Exf, bool bInt, bool bUpperLimit){
  // 'this' is a special pointer to the "current instance"
  this->index = index;
  this->E_cm = E_cm;
  this->dE_cm = dE_cm;
  this->wg = wg;
  this->dwg = dwg;
  this->Jr = Jr;
  for(int i=0; i<3; i++){
    this->G[i] = G[i];
    this->dG[i] = dG[i];
    this->L[i] = L[i];
    this->PT[i] = PT[i];
    this->dPT[i] = dPT[i];
  }
  this->Exf = Exf;
  this->bInt_flag = bInt;
  this->bUpperLimit = bUpperLimit;

  //cout << "made a resonance" << endl;
}

Resonance::~Resonance(){}

void Resonance::makeSamples(std::vector<std::vector<double> > Ref_sample, double smallestdE,
			    double smallestdwg, double smallestdG[3]){

  double mu, sigma;
  
  // First get the energy samples
  double corr = smallestdE/dE_cm;     // The correlation factor for this resonance energy
  E_sample.resize(NSamples);
  // Calculate correlated energies
  for(int s=0; s<NSamples; s++){
    double x2 = gsl_ran_gaussian(r, 1.0);
    x2 = corr*Ref_sample[s][0] + x2*sqrt(1. - gsl_pow_2(corr));
    E_sample[s] = E_cm + x2*dE_cm;
  }

  /*
  std::cout << "index: " << index << "\n";
  if(index==2 && bUpperLimit==false){
    for(int s=0; s<NSamples; s++)
      testfile << E_sample[s] << "\n";
  }
  */

  // ------------------------------
  // Now the resonance strength if it's known
  if(wg > 0.0){
    wg_sample.resize(NSamples);

    // Convert input values (expectation value and standard deviation)
    // to lognormal mu and sigma
    logNormalize(wg, dwg, mu, sigma);
    // Correlation parameter
    corr = smallestdwg*wg/dwg;

    for(int s=0; s<NSamples; s++){
      double x2 = gsl_ran_gaussian(r, 1.0);
      // Correlate everything with the first known resonance strength, wg_0
      // using: wg_i = p*wg_0 + wg_i*sqrt(1-p^2)
      x2 = corr*Ref_sample[s][1] + x2*sqrt(1.-gsl_pow_2(corr));
      // Convert normally distributed sample into lognormal for this resonance
      wg_sample[s] = gsl_sf_exp(mu + sigma*x2);
    }

    /*    
    std::cout << "index: " << index << "\n";
    if(index==2 && bUpperLimit==false){
      for(int s=0; s<NSamples; s++)
	testfile << wg_sample[s] << "\n";
    }
    */
    
  }

  //------------------------------
  // Now the partial widths
  for(int channel=0; channel<3; channel++){
    std::vector<double> G_temp;
    G_temp.resize(NSamples);
    // G_sample.resize(NSamples);
    // First if this channel is NOT an upper limit
    if(dG[channel] > 0.0){
      logNormalize(G[channel], dG[channel], mu, sigma);
      corr = smallestdG[channel] * G[channel]/dG[channel];
      for(int s=0; s<NSamples; s++){
	double x2 = gsl_ran_gaussian(r,1.0);
	x2 = corr*Ref_sample[s][channel+1] + x2*sqrt(1.-gsl_pow_2(corr));
	G_temp[s] = gsl_sf_exp(mu + x2*sigma);
      }

    } else {       // Or if it is an upper limit


    }
    
    G_sample.push_back(G_temp);
  } // end for(int channel=0; channel<3; channel++)

  /*
  if(index==0 && bUpperLimit==false){
    for(int s=0; s<NSamples; s++)
      testfile << G_sample[0][s] << "  " << G_sample[1][s] << "  " << G_sample[2][s] << "\n";
  }
  */

  
}


void Resonance::print(){

  //  cout << "--------------------------------------------------" << "\n";
  //cout << "     This is resonance: " << index << "\n";
  cout << " Resonace " << std::setw(3) << index <<"    E_cm = " << E_cm  << " +/- " << dE_cm << "\n";
  cout << "                 wg   = " << wg << " +/- " << dwg << "\n";
  cout << "                 Jr   = " << Jr << "\n";
  cout << "                 G1   = " << G[0] << " +/- " << dG[0] << " (L = " << L[0] << ")\n";
  if(bUpperLimit)
    cout << "                 PT   = " << PT[0] << " +/- " << dPT[0] << "\n";
  cout << "                 G2   = " << G[1] << " +/- " << dG[1] << " (L = " << L[1] << ")\n";
  if(bUpperLimit)
    cout << "                 PT   = " << PT[1] << " +/- " << dPT[1] << "\n";
  cout << "                 G3   = " << G[2] << " +/- " << dG[2] << " (L = " << L[2] << ")\n";
  if(bUpperLimit)
    cout << "                 PT   = " << PT[2] << " +/- " << dPT[2] << "\n";
  cout << "                 Exf  = " << Exf << "\n";
  cout << "           Integrated = " << bInt_flag << "\n";
  cout << "          Upper Limit = " << bUpperLimit << "\n";
  //  cout << "--------------------------------------------------" << "\n";
  cout << "\n";

}
void Resonance::write(){

  //  logfile << "--------------------------------------------------" << "\n";
  //logfile << "     This is resonance: " << index << "\n";
  logfile << " Resonace " << std::setw(3) << index <<"    E_cm = " << E_cm  << " +/- " << dE_cm << "\n";
  logfile << "                 wg   = " << wg << " +/- " << dwg << "\n";
  logfile << "                 Jr   = " << Jr << "\n";
  logfile << "                 G1   = " << G[0] << " +/- " << dG[0] << " (L = " << L[0] << ")\n";
  if(bUpperLimit)
    logfile << "                 PT   = " << PT[0] << " +/- " << dPT[0] << "\n";
  logfile << "                 G2   = " << G[1] << " +/- " << dG[1] << " (L = " << L[1] << ")\n";
  if(bUpperLimit)
    logfile << "                 PT   = " << PT[1] << " +/- " << dPT[1] << "\n";
  logfile << "                 G3   = " << G[2] << " +/- " << dG[2] << " (L = " << L[2] << ")\n";
  if(bUpperLimit)
    logfile << "                 PT   = " << PT[2] << " +/- " << dPT[2] << "\n";
  logfile << "                 Exf  = " << Exf << "\n";
  logfile << "           Integrated = " << bInt_flag << "\n";
  logfile << "          Upper Limit = " << bUpperLimit << "\n";
  //  logfile << "--------------------------------------------------" << "\n";
  logfile << "\n";

}
