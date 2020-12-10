#ifndef _Resonance_h_
#define _Resonance_h_

#include <iostream>

class Resonance{

 public:

  Resonance();
  ~Resonance();

  // Getters
  //  void getName();

  // Setters
  //void setName(std::string a){Name=a;}

 private:

  // Resonance control and bookkeeping
  bool bInt_flag, bisECorr;
  bool ErrorFlag;
  
  // Resonance parameters
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
