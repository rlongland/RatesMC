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
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <numeric>


#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_integration.h>

#include "Resonance.h"
#include "Reaction.h"
#include "Utilities.h"

using std::cout;
using std::endl;

Reaction::Reaction(){

  smallestdE = 0.0;
  smallestdwg = 0.0;
  for(int i=0; i<3; i++)
    smallestdG[i] = 0.0;

}

Reaction::~Reaction(){}

void Reaction::setNonResonant(double s, double sp, double spp, double ds, double cutoffe, int part){
  S[part] = s;
  Sp[part] = sp;
  Spp[part] = spp;
  dS[part] = ds;
  CutoffE[part] = cutoffe;

  std::vector<double> A(NSamples);
  ARate.push_back(A);
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
  Resonance Res(*this, i, E_cm, dE_cm, wg, dwg, Jr,
		G, dG, L, PT, dPT, Exf, bInt, bUpperLimit);
  //Res.print();
  // Add that resonance to the list of resonances
  Resonances.push_back(Res);
}

void Reaction::printName(){
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
  cout << "--------------------------------------------------" << "\n";
  

}

void Reaction::writeReaction(){

  logfile << "--------------------------------------------------" << "\n";
  logfile << "     This is reaction: " << Name << "\n";
  logfile << "   Z0        Z1       Z2" << "\n";
  logfile << std::setw(5) << Z0 << "     " << std::setw(5) << Z1 << "    " << std::setw(5) << Z2 << "\n";
  logfile << "   M0        M1       M2" << "\n";
  logfile << M0 << "   " << M1 << "   " << M2 << "\n";
  logfile << "   J0        J1       J2" << "\n";
  logfile << J0 << "   " << J1 << "   " << J2 << "\n";
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
    //    std::cout << "Res: " << res->getIndex() << "\n";
    res->write();
  }
  //std::cout << "Done!\n";

}

//----------------------------------------------------------------------
// Set up the contribution file header
void Reaction::setupContribHeader(){

  contribfile << "#Temp" << "   ";
  for(int i=0;i<3;i++)
    contribfile << "A-Rate-1   ";
  for(int i=0;i<3;i++)
    contribfile << "A-Rate-2   ";
  for(Resonance Res : Resonances){
    for(int i=0; i<3; i++){
      contribfile << "Res" << Res.getIndex()+1 << "       ";
    }
  }
  contribfile << endl;

}

//----------------------------------------------------------------------
// Calculate the resonant reaction rate
double Reaction::calcResonant(double Temp){

  // Factor for final rate (iliadis Eqn 3.114)
  double mue = M0*M1/(M0+M1);
  double RateFactorNarrow = (1.5399e11/pow(mue*Temp,1.5));
  double RateFactorBroad = 3.7318e10/(pow(mue,0.5)*pow(Temp,1.5));
  
  // Vector of rate storage
  std::vector<double> individualRate;   // individualRate stores the "classical" rate for each resonance
  double classicalRate=0.0;             // classicalRate is to summed individialRate

  // Loop over all resonances
  for(Resonance &R : Resonances){
    // if the resonance is narrow
    if(!R.getisBroad()){
      //      std::cout << "Resonance " << R.getIndex() << " at "
      //      		<< R.getE_cm() << " keV is narrow\n";
      individualRate.push_back(RateFactorNarrow*R.calcNarrow(Temp));
      // Multiply every sample by the reaction rate constant above
      R.scaleByConstant(RateFactorNarrow);
      // If it's broad
    } else {
      //      std::cout << "Resonance " << R.getIndex() << " at "
      //      		<< R.getE_cm() << " keV is being numerically integrated\n";
      individualRate.push_back(RateFactorBroad*R.calcBroad(Temp));
      // Multiply every sample by the reaction rate constant above
      R.scaleByConstant(RateFactorBroad);
    }


    
    //    std::transform(individialRate.begin(), individialRate.end(), individialRate.begin(),
    //		   [RateFactor](int &c){ return c*RateFactor; });

    // The classical rate is the last element of the individual rates
    //    classicalRate += individialRate.back();

    // Once calculated, print the rate and some diagnostics for each resonance
    R.printRate();
  }

  classicalRate = std::accumulate(individualRate.begin(), individualRate.end(), 0.0);

  std::cout << "Total classical rate from resonances = " << classicalRate << "\n";

  return classicalRate;
}

//----------------------------------------------------------------------
// Collect the resonant rates for sample s
std::vector<double> Reaction::getResonantRateSample(int s){

  std::vector<double> Rate_s;
  
  for(Resonance R : Resonances){
    //R.printRate();

    Rate_s.push_back(R.getRateSample(s));
  }
  //std::cout << Rate_s[0] << "   ";
  return Rate_s;
}

//----------------------------------------------------------------------
// Calculate the direct capture, non-resonant part of the reaction rate
double Reaction::calcNonResonant(double Temp, int j){

  //  std::cout << Temp << " " << j << "\n";
  double mue = M0*M1/(M0+M1);
  //  double mu, sigma;
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
  

  /*cout << C1e << "\t" << C2e << "\t" << C3e << "\t" << C4e << "\t" <<
    C5e << "\t" << C6e << "\t" << C7e << endl;
    cout << "Cutoff T = " << cutoff_T << endl;*/

  // Calculate the rate first and then sample at the end.
  
  // Turn off the error handler in case exp returns zero
  gsl_set_error_handler_off();
  ADRate = 1.0e-3*(C1e/pow(Temp,2./3.))*gsl_sf_exp(-C2e/pow(Temp,1./3.)-
					    pow(Temp/cutoff_T,2.))*
    (1 + C3e*pow(Temp,1./3.) + C4e*pow(Temp,2./3.) + C5e*Temp +
     C6e*pow(Temp,4./3.) + C7e*pow(Temp,5./3.));

  gsl_set_error_handler(NULL); // Turn error handler back on again
  
  
  /*cout << (C1e/pow(Temp,2./3.))*exp(-C2e/pow(Temp,1./3.)-
    pow(Temp/cutoff_T,2.)) << endl;
    cout << C3e*pow(Temp,1./3.) << endl;
    cout << C4e*pow(Temp,2./3.) << endl;
    cout <<  C5e*Temp << endl;

    cout << "ADRate = " << ADRate << endl;
    cout << "About to Lognormalise" << endl;*/

  // If ADRate is negative, set to zero, this is unphysical
  if(ADRate < 0.){
    ErrorFlag = true;
    logfile << "\tWARNING: The non-resonant part caused a negative rate, \n\t\tsetting to zero." << endl;
    ADRate = 0.0;
    for(int i=0;i<NSamples;i++){
      ARate[j][i]=0.0;
    }
  } else {
    logNormalize(ADRate,dS[j]*ADRate,mu,sigma);
    
    // if mu and sigma return as zero, the rate is zero (very small)
    if((mu==0. && sigma==0.)){
      ADRate=0.0;
      for(int i=0;i<NSamples;i++){
	ARate[j][i]=0.0;
      }
    } else {
      for(int i=0;i<NSamples;i++){
	ARate[j][i] = gsl_ran_lognormal(r,mu,sigma);///(1.5399e11/pow(mue*Temp,1.5));
	//	std::cout << ARate[j][i] << "\n";
      }
    }
  }


  ADRate = ADRate;///(1.5399e11/pow(mue*Temp,1.5));
  
  return ADRate;
}

//----------------------------------------------------------------------
// Calculate the non-resonant part of the reaction rate by integrating the
// astrophysical s-factor
// Ugly-ass hack
Reaction * ReactionPtr;
double NonResonantIntegrandWrapper(double x, void *params)
{
  return ReactionPtr->NonResonantIntegrand(x, params);
}
// end ugly-ass hack
double Reaction::calcNonResonantIntegrated(double Temp, int j){

  //  std::cout << Temp << " " << j << "\n";
  double mue = M0*M1/(M0+M1);
  //  double mu, sigma;

  double ADRate,mu,sigma;

  // Calculate the rate first and then sample at the end.
  double E_min = 0.0;
  double E_max = CutoffE[j]/1000.0;
  //  std::cout << E_max << std::endl;
    // GSL Integration functions
  double result, error;
  
  // Define integration limits
  double x = E_min, x1 = E_max;

  // The array that needs to be passed to the integration function
  double alpha[7];
  alpha[0] = j;  // The index of the non-resonant part
  alpha[1] = mue;
  alpha[2] = Temp;
  alpha[3] = S[j]/1000.0;
  alpha[4] = Sp[j];
  alpha[5] = Spp[j]*1000.0;
  alpha[6] = Z0*Z1;
  //alpha[1] = (int)writeIntegrand;


  // Turn off the error handler
  gsl_set_error_handler_off();

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_function F;
  // Can't use Integrand directly because GSL is shit
  F.function = &NonResonantIntegrandWrapper;
  //  F.function = &Integrand;
  F.params = &alpha;

  int status = gsl_integration_qag(&F,      // Function to be integrated
				   E_min,   // Start of integration
				   E_max,   // End of integration
				   0,       // absolute error
				   1e-3,    // relative error
				   1000,    // max number of steps (cannot exceed size of workspace
				   6,       // key - (6=61 point Gauss-Kronrod rules)
				   w,       // workspace
				   &result, // The result
				   &error);


  gsl_set_error_handler(NULL);

  ADRate = result*(3.7318e10/(sqrt(mue)*pow(Temp,1.5)));

  //std::cout << "AD Rate = " << ADRate << std::endl;
  

  // If ADRate is negative, set to zero, this is unphysical
  if(ADRate < 0.){
    ErrorFlag = true;
    logfile << "\tWARNING: The non-resonant part caused a negative rate, \n\t\tsetting to zero." << endl;
    ADRate = 0.0;
    for(int i=0;i<NSamples;i++){
      ARate[j][i]=0.0;
    }
  } else {
    // Either use the fractional uncertainty in S
    if(dS[j] >= 0){
      logNormalize(ADRate,dS[j]*ADRate,mu,sigma);
      // Or interpret it as a factor uncertainty
    } else {
      mu = gsl_sf_log(ADRate);
      sigma = gsl_sf_log(-dS[j]);
    }
    // if mu and sigma return as zero, the rate is zero (very small)
    if((mu==0. && sigma==0.)){
      ADRate=0.0;
      for(int i=0;i<NSamples;i++){
	ARate[j][i]=0.0;
      }
    } else {
      for(int i=0;i<NSamples;i++){
	ARate[j][i] = gsl_ran_lognormal(r,mu,sigma);///(1.5399e11/pow(mue*Temp,1.5));
	//	std::cout << ARate[j][i] << "\n";
      }
    }
  }


  ADRate = ADRate;///(1.5399e11/pow(mue*Temp,1.5));
  
  return ADRate;
}

double Reaction::NonResonantIntegrand(double x, void * params){

  //  std::cout << "inside " << x << std::endl; 
  double* par = (double*)params;
  int index = (int)par[0];
  double mue = (double)par[1];
  double T = (double)par[2];
  double S = (double)par[3];
  double Sp = (double)par[4];
  double Spp = (double)par[5];
  double Z0Z1 = (double)par[6];

  //  std::cout << index << " " << mue << " " << T << std::endl;
  //  std::cout << S << std::endl;
  double Ssum = S + Sp*x + 0.5*Spp*x*x;

  //  std::cout << "Ssum = " << Ssum << std::endl;
  
  double eta = 0.989510*Z0Z1*sqrt(mue/x);
  double Sommerfeld = gsl_sf_exp(-eta);
  double Boltzmann = gsl_sf_exp(-11.605*x/T);
  
  double integrand = Ssum*Sommerfeld*Boltzmann;

  return integrand;
}

//----------------------------------------------------------------------
// Prepare the Monte Carlo samples. 
void Reaction::prepareSamples(){

  cout << "Preparing " << NSamples << " samples\n\n";

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
  for(Resonance &Res : Resonances){
    //std::cout << Res.getIndex() << "\n";
    Res.makeSamples(Ref_sample, smallestdE, smallestdwg, smallestdG);
    //std::cout << "done\n";
  }

  
}


//----------------------------------------------------------------------
// Write all input parameter samples to a file
void Reaction::writeSamples(){

  std::cout << "Writing sample file. This may take a while...\n";
  
  std::ofstream samplefile;
  samplefile.open("ParameterSamples.dat");
  
  // Each column corresponds to a parameter, rows are samples
  // Write the header
  samplefile <<    "                                           ";
  for(Resonance Res : Resonances){
    samplefile << " |         Resonance " << std::setw(3) << Res.getIndex()
	       << " at E_cm = " << std::setw(7) << Res.getE_cm() << " MeV          ";
  }
  samplefile << "\n";
  //               "1234567890x1234567890x1234567890x123456789012"
  samplefile <<    " Standard1  Standard2  Standard3  Standard4";
  for(Resonance Res : Resonances){
    samplefile << " |         E         wg         G1         G2         G3";
  }
  samplefile << "\n";
  
  // Now the samples
  for(int s=0; s<NSamples; s++){
    //    samplefile << "  ";
    // Reference randoms for correlations
    for(int channel=0; channel<4; channel++){
      samplefile << std::setw(10) << Ref_sample[s][channel] << " ";
    }
    samplefile << " ";
    // Now the resonances
    for(Resonance Res: Resonances){
      Res.writeSamples(samplefile, s);
    }
    
    samplefile << "\n";
  }

  samplefile.close();
  
  std::cout << "Done!\n\n";
}
