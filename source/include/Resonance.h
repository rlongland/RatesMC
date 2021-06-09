#ifndef _Resonance_h_
#define _Resonance_h_

#include <iostream>

class Resonance{

 public:

  // Constructor
  Resonance(int index=0, double E_cm=0.0, double dE_cm=0.0, double wg=0.0, double dwg=0.0, double Jr=0.0,
	    double G[3]={}, double dG[3]={}, int L[3]={0}, double PT[3]={}, double dPT[3]={},
	    double Exf=0, bool bInt=false, bool bUpperLimit=false); 
  // Destructor
  ~Resonance();

  // Getters
  void getIndex();
  void getE_cm();

  // Setters
  void setIndex(int i){index=i;}

  // print a summary of the resonance
  void print();
  void write();
  
 private:

  // Resonance control and bookkeeping
  bool bInt_flag, bisECorr, bUpperLimit;
  bool ErrorFlag;
  
  // Resonance parameters
  int index;
  int L[3];
  double E_cm, dE_cm, wg, dwg, Jr, G[3], dG[3], Exf, PT[3], dPT[3];
  bool bwg_known;
  bool bUseInrate;
  
  // Correlations
  int CorresRes;
  double Ref_sample[3], ERef_sample, Frac;
  double correlation, Gcorrelations, Ecorrelation;
};




#endif
