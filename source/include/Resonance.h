#ifndef _Resonance_h_
#define _Resonance_h_

#include <iostream>
#include <vector>

class Reaction;

class Resonance {

 public:

  // Constructor
  Resonance(Reaction & R,
	    int index=0, double E_cm=0.0, double dE_cm=0.0, double wg=0.0, double dwg=0.0, double Jr=0.0,
	    double G[3]={}, double dG[3]={}, int L[3]={0}, double PT[3]={}, double dPT[3]={},
	    double Exf=0, bool bInt=false, bool bUpperLimit=false); 
  // Destructor
  ~Resonance();

  // Getters
  int getIndex(){return index;}
  void getE_cm();

  // Setters
  void setIndex(int i){index=i;}

  void makeSamples(std::vector<std::vector<double> > Ref_sample, double smallestdE,
		   double smallestdwg, double smallestdG[3]);
  void writeSamples(std::ofstream& samplefile, int s);
  
  // print a summary of the resonance
  void print();
  void write();
  
 private:

  // Access to the reaction
  Reaction & Reac;
  
  // Resonance control and bookkeeping
  bool bInt_flag, bisECorr, bUpperLimit;
  bool ErrorFlag;
  
  // Resonance parameters
  int index;
  int L[3];
  double E_cm, dE_cm, wg, dwg, Jr, G[3], dG[3], Exf, PT[3], dPT[3];
  double M0,M1,M2,J0,J1,J2;
  int Z0,Z1,Z2;
  
  //  bool bwg_known;
  bool bUseInrate;

  // Sampled parameters
  std::vector<double> E_sample;
  std::vector<double> wg_sample;
  std::vector<std::vector<double> > G_sample;
  
  // Correlations
  int CorresRes;
  double Ref_sample[3], ERef_sample, Frac;
  double correlation, Gcorrelations, Ecorrelation;
};




#endif
