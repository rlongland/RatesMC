#ifndef _Reaction_h_
#define _Reaction_h_

#include <iostream>
#include <vector>
#include "Resonance.h"

class Reaction{

 public:

  Reaction();
  ~Reaction();

  // Getters
  void getName();

  // Setters
  void setName(std::string a){Name=a;}
  void setMasses(double m0, double m1, double m2){M0=m0; M1=m1; M2=m2;}
  void setCharges(int z0, int z1, int z2){Z0=z0; Z1=z1; Z2=z2;}
  void setSpins(int j0, int j1, int j2){J0=j0; J1=j1; J2=j2;}
  void setSeparationEnergies(double qin, double qout){Q=qin; Qexit=qout;}
  void setR0(double r0){R0=r0;}
  void setGammaIndex(int gindex){Gamma_index = gindex;}
  void setNonResonant(double, double, double, double, double, int);
  void addResonance(int, double, double, double, double, double,
		    double, double, int, double, double,
		    double, double, int, double, double,
		    double, double, int, double, double,
		    double, bool, bool);
  double calcNonResonant(double Temp, int j);
  void prepareSamples();

  // Getters
  double getsmallestdE(){return smallestdE;}
  int getGamma_index(){return Gamma_index;}
  
  // Print a summary of the reaction
  void printReaction();
  void writeReaction();

  double M0,M1,M2,J0,J1,J2,Q,dQ,Qexit,dQexit,R0;
  int Z0,Z1,Z2;
  
 private:

  std::string Name;
  double S[2],Sp[2],Spp[2],dS[2],CutoffE[2];
  int NRes,ULNRes,Gamma_index;
  int *NChannels;
  double smallestdE;
  double smallestdwg;
  double smallestdG[3];
  std::vector<Resonance> Resonances;

  // Reaction-wide Monte Carlo
  std::vector<std::vector<double> > Ref_sample;
  std::vector<std::vector<double> > DRate;
  
  
};



#endif
