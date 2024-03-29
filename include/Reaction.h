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

#ifndef _Reaction_h_
#define _Reaction_h_

#include <iostream>
#include <vector>
#include <array>
#include "Resonance.h"
#include "Interference.h"

class Reaction{

 public:

  Reaction();
  ~Reaction();

  // Getters
  void printName();

  // Setters
  void setName(std::string a){Name=a;}
  void setMasses(double m0, double m1, double m2){M0=m0; M1=m1; M2=m2;}
  void setCharges(int z0, int z1, int z2){Z0=z0; Z1=z1; Z2=z2;}
  void setSpins(double j0, double j1, double j2){J0=j0; J1=j1; J2=j2;}
  void setSeparationEnergies(double qin, double qout){Q=qin; Qexit=qout;}
  void setR0(double r0){R0=r0;}
  void setGammaIndex(int gindex){Gamma_index = gindex;}
  void setNonResonant(double, double, double, double, double, int);
  void setNonResonantTable(std::vector<double>, std::vector<double>,
													 std::vector<double>, int);
	void setCutoffE(double EMax, int part){CutoffE[part] = EMax;}
  void addResonance(int, double, double, double, double, double,
										double, double, int, double, double,
										double, double, int, double, double,
										double, double, int, double, double,
										double, bool, bool, bool, bool, int, double);
  void addInterference(int,int);
	void addResonanceToInterference(int,
																	double, double,	double,
																	double, double, int, double, double,
																	double, double, int, double, double,
																	double, double, int, double, double,
																	double, bool, int);
  double calcResonant(double Temp);
  double calcNonResonant(double Temp, int j);
  double calcNonResonantIntegrated(double Temp, int j);
  double calcNonResonantTabulated(double Temp, int j);
  double NonResonantIntegrand(double x, void * params);
  double NonResonantTabulatedIntegrand(double x, void * params);
  void prepareSamples();
  void writeSamples();
  void writeSFactor(bool MCSamples);
  void setupSFactorHeader(std::ofstream &sfactorfile);

	void CombineResonancePossibilities();
	
  // Getters
  std::string getName(){return Name;}
  double getsmallestdE(){return smallestdE;}
  int getGamma_index(){return Gamma_index;}
	Resonance getLastResonance(){return Resonances.back();}
	int getNResonances(){return Resonances.size();}

  // Get rates
  // Analytical Rate
  double getARate(int s, int index){return ARate[index][s];}
  // Resonant rate (vector's length is No. of resonances)
  std::vector<double> getResonantRateSample(int s);

  // Set up the contribution file header
  void setupContribHeader();
  
  // Print a summary of the reaction
  void printReaction();
  void writeReaction();

  double M0,M1,M2,J0,J1,J2,Q,dQ,Qexit,dQexit,R0;
  int Z0,Z1,Z2;

	double smallestdE, EofSmallestdE;
  double smallestdwg, EofSmallestdwg;
  double smallestdG[3], EofSmallestdG[3];

 private:

  std::string Name;
  double S[2],Sp[2],Spp[2],dS[2],CutoffE[2];
  int Gamma_index;

	// Tabulated S-factor
	// std::vector<std::array<double, 2> > SFactorE;
	// std::vector<std::array<double, 2> > SFactorS;
	// std::vector<std::array<double, 2> > SFactordS;
	std::array<std::vector<double>, 2> SFactorE;
	std::array<std::vector<double>, 2> SFactorS;
	std::array<std::vector<double>, 2> SFactordS;

  std::vector<Resonance> Resonances;
	std::vector<Interference> Interferences;

  // Reaction-wide Monte Carlo
  std::vector<std::vector<double> > Ref_sample;
	std::vector<std::vector<double> > ARate;
	std::vector<std::vector<double> > SScale_sample;
  
};



#endif
