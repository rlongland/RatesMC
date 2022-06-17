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

#ifndef _Interference_h_
#define _Interference_h_

#include <iostream>
#include <vector>


class Reaction;
class Resonance;

class Interference {

 public:


  // Constructor
  Interference(Reaction & R, int i, int IntfSign); 
  // Destructor
  ~Interference();

  // Getters
  int getIndex(){return index;}
  
  // Setters
  void setIndex(int i){index=i;}

	void addResonance(Resonance* Res, int index);

	// Make the MC samples
	void makeSamples(std::vector<std::vector<double> > Ref_sample,
									 double, double, double[3]);
	
  // Functions that do stuff
  // Functions to calculate the rate from this resonance. They return
  // the traditional rate and fill 'Rate', which is a vector of rate
  // samples for this resonance at temperature T
  double calcBroad(double T);

  void scaleByConstant(double RateFactor);
  
  void printRate();

	void putRateSample(int s, double rate){Rate_sample[s] = rate;}
  double getRateSample(int s){return Rate_sample[s];}

	double getSFactor(double E);
	
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
  
  // Control and bookkeeping
  bool ErrorFlag;

	// TODO Should just hold two resonances, not the rest of this junk
	Resonance* Res[2];

	int IntfSign;
	
	// Resonance parameters
  int index;
  int L[3];
  double E_cm, dE_cm, Jr, G[3], dG[3], Exf, PT[3], dPT[3];
  double M0,M1,M2,J0,J1,J2;
  int Z0,Z1,Z2;
  int NChannels;
  
  double R, mue;
  
  // Sampled parameters
  std::vector<double> E_sample;               // Energies
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
