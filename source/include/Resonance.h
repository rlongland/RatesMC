#ifndef _Resonance_h_
#define _Resonance_h_

#include <iostream>

class Resonance{

 public:

  // Constructor
  Resonance(int index=0, double E_cm=0.0, double dE_cm=0.0, double wg=0.0, double dwg=0.0); 
  // Destructor
  ~Resonance();

  // Getters
  void getIndex();
  void getE_cm();

  // Setters
  void setIndex(int i){index=i;}

  // print a summary of the resonance
  void print();
  
 private:

  // Resonance control and bookkeeping
  bool bInt_flag, bisECorr;
  bool ErrorFlag;
  
  // Resonance parameters
  int index;
  int L[3];
  double E_cm, dE_cm, wg, dwg, J, G, dG, Exf, PTMean[3], dPTMean[3];
  bool bwg_known;
  bool bUseInrate;
  
  // Correlations
  int CorresRes;
  double Ref_sample[3], ERef_sample, Frac;
  double correlation, Gcorrelations, Ecorrelation;
};




#endif
