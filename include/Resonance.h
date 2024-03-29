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

#ifndef _Resonance_h_
#define _Resonance_h_

#include <iostream>
#include <vector>


class Reaction;

static double defaultdinit[3] = {0.0,0.0,0.0};
static int defaultiinit[3] = {0,0,0};

class Resonance {

 public:


  // Constructor
  Resonance(Reaction & R,
	    int index=0, double E_cm=0.0, double dE_cm=0.0, double wg=0.0, double dwg=0.0, double Jr=0.0,
						double G[3]=defaultdinit, double dG[3]=defaultdinit, int L[3]=defaultiinit,
						double PT[3]=defaultdinit, double dPT[3]=defaultdinit,
						double Exf=0, bool bInt=false, bool bUpperLimit=false,
						bool isECorrelated=false, bool isWidthCorrelated=false,
						int CorresRes=0, double Frac=1.0); 
  // Destructor
  ~Resonance();

  // Getters
  int getIndex(){return index;}
  double getE_cm(){return E_cm;}
  double getdE_cm(){return dE_cm;}
  double getisBroad(){return isBroad;}
	int getCorresRes(){return CorresRes;}
	double getFrac(){return Frac;}
	double getExf(){return Exf;}
	double getESample(int s){return E_sample[s];}
	double getGSample(int i,int s){return G_sample[i][s];}
	double geterFrac(int i,int s){return erFrac[i][s];}
	double getG(int i){return G[i];}
	int getL(int i){return L[i];}
	double getM0(){return M0;}
	double getM1(){return M1;}
	double getM2(){return M2;}
	int getZ0(){return Z0;}
	int getZ1(){return Z1;}
	int getZ2(){return Z2;}
	double getR(){return R;}
	double getmue(){return mue;}
	int getNChannels(){return NChannels;}
	double getJr(){return Jr;}
	
	// Setters
  void setIndex(int i){index=i;}

  // Functions that do stuff
  void makeSamples(std::vector<std::vector<double> > Ref_sample, double smallestdE,
		   double smallestdwg, double smallestdG[3]);
  void writeSamples(std::ofstream& samplefile, int s);
  // Functions to calculate the rate from this resonance. They return
  // the traditional rate and fill 'Rate', which is a vector of rate
  // samples for this resonance at temperature T
  double calcBroad(double T);
  double calcNarrow(double T);

  void scaleByConstant(double RateFactor);
  
  void printRate();

	void putRateSample(int s, double rate){Rate_sample[s] = rate;}
  double getRateSample(int s){return Rate_sample[s];}

	double getSFactor(double E, int samp);
	
  // Rate for a single narrow resonance
  double singleNarrow(double wg, double E, double T);
  // Function to integrate a broad resonance
  double NumericalRate(double T,
		       double E, double G0, double G1, double G2,
		       double erFrac0, double erFrac1, double erFrac2,
		       bool writeIntegrand);
  // The resonance integrand to be integrated in the function above
  //extern "C" {
  double Integrand(double x, void * params);
	int rhs(double x, const double y[], double dydx[], void * params);
  //}
  
  // print a summary of the resonance
  void print();
  void write();
  
 private:

  // Access to the reaction
  Reaction & Reac;
  
  // Resonance control and bookkeeping
  bool isBroad, isECorrelated, isWidthCorrelated, isUpperLimit;
  bool ErrorFlag;
  
  // Resonance parameters
  int index;
  int L[3];
  double E_cm, dE_cm, wg, dwg, Jr, G[3], dG[3], Exf, PT[3], dPT[3];
  double M0,M1,M2,J0,J1,J2;
  int Z0,Z1,Z2;
  int NChannels;
	int CorresRes;
	double Frac;
  
  double R, R_exit, mue;
  
  //  bool bwg_known;
  //bool bUseInrate;

  // Sampled parameters
  std::vector<double> E_sample;               // Energies
  std::vector<double> wg_sample;              // Strengths
  std::vector<std::vector<double> > G_sample; // Partial widths
  std::vector<std::vector<double> > erFrac;   // Energy shift effect on widths
  
  // The rate from this resonance
  double classicalRate;
  std::vector<double> Rate_sample;   // length is number of samples
  double MeanRate, MedianRate, RateMu, RateSigma;

};


// Junk to get GSL to play with classes
/*
template< typename F >
class gsl_function_pp : public gsl_function {
public:
  gsl_function_pp(const F& func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params=this;
  }
private:
  const F& _func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->_func(x);
  }
};
*/

#endif
