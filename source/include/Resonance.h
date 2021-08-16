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
  double getE_cm(){return E_cm;}
  double getisBroad(){return isBroad;}
  
  // Setters
  void setIndex(int i){index=i;}

  // Functions that do stuff
  void makeSamples(std::vector<std::vector<double> > Ref_sample, double smallestdE,
		   double smallestdwg, double smallestdG[3]);
  void writeSamples(std::ofstream& samplefile, int s);
  // Functions to calculate the rate from this resonance. They return
  // the traditional rate and fill 'Rate', which is a vector of rate
  // samples for this resonance at temperature T
  double calcBroad(double T, std::vector<double> &Rate);
  double calcNarrow(double T, std::vector<double> &Rate);

  double singleNarrow(double wg, double E, double T);
  
  // print a summary of the resonance
  void print();
  void write();
  
 private:

  // Access to the reaction
  Reaction & Reac;
  
  // Resonance control and bookkeeping
  bool isBroad, isECorr, isUpperLimit;
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
  std::vector<double> E_sample;               // Energies
  std::vector<double> wg_sample;              // Strengths
  std::vector<std::vector<double> > G_sample; // Partial widths
  std::vector<std::vector<double> > erFrac;   // Energy shift effect on widths
  
  // Correlations
  int CorresRes;
  double Ref_sample[3], ERef_sample, Frac;
  double correlation, Gcorrelations, Ecorrelation;
};




#endif
