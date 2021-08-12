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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>

#include "Resonance.h"
#include "Utilities.h"

using std::cout;
using std::endl;

Resonance::Resonance(Reaction & R,
		     int index, double E_cm, double dE_cm, double wg, double dwg, double Jr,
		     double G[3], double dG[3], int L[3], double PT[3], double dPT[3],
		     double Exf, bool bInt, bool bUpperLimit) : Reac(R){
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

  this->M0 = R.M0;
  this->M1 = R.M1;
  this->M2 = R.M2;
  this->Z0 = R.Z0;
  this->Z1 = R.Z1;
  this->Z2 = R.Z2;
  
  //  std::cout << "The Gamma_index is " << Reac.getGamma_index() << "\n";

  //cout << "made a resonance" << endl;
}

Resonance::~Resonance(){}

void Resonance::makeSamples(std::vector<std::vector<double> > Ref_sample, double smallestdE,
			    double smallestdwg, double smallestdG[3]){

  double mu, sigma;
  double P;

  double mue = M0*M1/(M0+M1);
  double R = Reac.R0*(pow(M0,(1./3.))+pow(M1,(1./3.)));

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

    // check that the uncertainty is reasonable
    if(isZero(dwg)){
      std::cout << "ERROR: You MUST specify a resonance strength uncertainty for resonance: " <<
	index << " at " << E_cm << " keV\n";
      std::cout << "       wg = " << wg << " +/- " << dwg << "\n";
      std::exit(EXIT_FAILURE);
    }
    
    wg_sample.resize(NSamples);

    // Convert input values (expectation value and standard deviation)
    // to lognormal mu and sigma
    logNormalize(wg, dwg, mu, sigma);

    // Correlation parameter
    corr = smallestdwg*wg/dwg;

    // Generate the correlated samples
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
    
  } else {

    //------------------------------
    // Now the partial widths
    for(int channel=0; channel<3; channel++){
      // Skip everything if G[channel]=0
       
      // G_sample.resize(NSamples);
      std::vector<double> G_temp;
      G_temp.resize(NSamples);

      // Skip channel 2 (third channel) if it's zero, otherwise make
      // sure the partial width is defined for G1 and G2
      if(channel == 2 && isZero(G[channel])){
	continue;
      } else if(channel < 2 && isZero(G[channel]) ){
	std::cout << "ERROR: You MUST specify a partial width for resonance: " <<
	  index << " at " << E_cm << " keV\n";
	std::cout << "       G" << channel+1 << " = " << G[channel] << " +/- " << dG[channel] << "\n";
	std::exit(EXIT_FAILURE);
      }  

      // First if this width is a normal, known width
      if(!(bUpperLimit && isZero(dG[channel]))){
	
	// check that the uncertainty is reasonable
	if(isZero(dG[channel])){
	  std::cout << "ERROR: You MUST specify partial width uncertainties for resonance: " <<
	    index << " at " << E_cm << " keV\n";
	  std::cout << "       G" << channel+1 << " = " << G[channel] << " +/- " << dG[channel] << "\n";
	  std::exit(EXIT_FAILURE);
	}

	// Find the lognormal parameters to generate the random partial width
	logNormalize(G[channel], dG[channel], mu, sigma);
	
	// Calculate the correlated partial widths for this channel
	corr = smallestdG[channel] * G[channel]/dG[channel];
	for(int s=0; s<NSamples; s++){
	  double x2 = gsl_ran_gaussian(r,1.0);
	  x2 = corr*Ref_sample[s][channel+1] + x2*sqrt(1.-gsl_pow_2(corr));
	  G_temp[s] = gsl_sf_exp(mu + x2*sigma);
	}

	// Or if it is an upper limit
      } else {     

	double A;  // constants

	// Is this the Gamma channel?
	if(channel == Reac.getGamma_index()){
	  std::cout << "Res " << index << " at E_cm = " << E_cm <<  " keV:\n" 
		    << "      Channel " << channel << " is Gamma_gamma\n";
	  A = (8.0*M_PI*(L[channel]+1)/(L[channel] *
					gsl_pow_2(gsl_sf_doublefact(2*L[channel]+1)))) *
	    gsl_pow_int((E_cm/197.326968),(2*L[channel]+1));
	  for(int s=0;s<NSamples;s++){
	    do{
	      // Randomise Porter Thomas according to it's uncertainty.
	      double PTi=0.0;
	      if(dPT[channel] < 0.0){
		// If it's negative, treat it as a factor uncertainty
		double mu = gsl_sf_log(PT[channel]);
		double sigma = gsl_sf_log(-dPT[channel]);
		PTi = gsl_ran_lognormal(r,mu,sigma);
	      } else {
		// otherwise just normal randomization
		double dof = 2.*gsl_pow_2(PT[channel]/dPT[channel]);
		PTi = PT[channel]*gsl_ran_chisq(r,dof)/dof;
		if(gsl_isinf(dof))PTi=PT[channel];
	      }
		
	      G_temp[s] = A*PTi*gsl_ran_chisq(r,1.0);
		
	      ptfile << index << "  " << s << "  " << PTi << "  " << G_temp[s]/A << endl;
		
	    } while(G_temp[s] > G[channel] && G[channel] != 0.0);

	  
	  }
	  // if it's not the Gamma channel
	} else { 
	  
	  std::cout << "Res " << index << " at E_cm = " << E_cm <<  " keV:\n"
		    << "      Channel " << channel << " is Gamma_particle\n";

	  // If it's an exit particle, we need to account for the Q-value
	  if(E_cm > 0.0 || channel != 0){
	    if(channel == 0){
	      P = PenFactor(E_cm,L[channel],M0,M1,Z0,Z1,R);
	      cout << E_cm << " " << L[channel] << " " << M0 << " " << M1 << " " << M2 << "\n";
	    } else if(channel == 1){
	      P = PenFactor(E_cm+Reac.Q-Reac.Qexit-Exf,L[channel],M0+M1-M2,M2,
			    Z0+Z1-Z2,Z2,R);
	    } else if(channel == 2){
	      P = PenFactor(E_cm+Reac.Q-Reac.Qexit,L[channel],M0+M1-M2,M2,
			    Z0+Z1-Z2,Z2,R);
	    }
	    // ((hbar c)^2)/u = 41.80159024 MeV^2 fm^2
	    A = 2.0*41.80161396*P/(mue*gsl_pow_2(R));
	  } else {
	    P = 1.0;
	    A = 1.0;
	  }
	  for(int s=0;s<NSamples;s++){
	    // Randomise Porter Thomas mean value
	    double PTi=0.0;
	    // Either as a factor uncertainty
	    if(dPT[channel] < 0.0){
	      double mu = gsl_sf_log(PT[channel]);
	      double sigma = gsl_sf_log(-dPT[channel]);
	      PTi = gsl_ran_lognormal(r,mu,sigma);
	      // or as a normal Gaussian uncertainty
	    } else {
	      double dof = 2.*gsl_pow_2(PT[channel]/dPT[channel]);
	      PTi = PT[channel]*gsl_ran_chisq(r,dof)/dof;
	      if(gsl_isinf(dof))PTi=PT[channel];
	    }
	    
	    // if an upper limit is entered, make sure the sample is
	    // less than it.
	    do{
	      G_temp[s] = A*PTi*gsl_ran_chisq(r,1.0);
	      
	      ptfile << s << "  " << index << "  " << P << "  " << PTi << "  " << G_temp[s] << endl;
	    } while(G_temp[s] > G[channel]);
	  } // for(int s=0; s<NSamples; s++)
	} // if(Gamma_gamma) else { 
      } // End else if it's an upper limit

      // By this point, we should have a bunch of samples for this
      // channel. Put them in the G_sample vector for saving
      G_sample.push_back(G_temp);
    } // end for(int channel=0; channel<3; channel++)
  } // end if(wg > 0) else
  
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
