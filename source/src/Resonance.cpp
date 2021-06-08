/* ======================================================================
   Resonance.cpp
   Author: R. Longland
   Date: 2021-01-27
   
   Description: Contains all of the resonance specific stuff
   ======================================================================
*/
#include <iostream>
#include <iomanip>
#include "stdio.h"
#include "Resonance.h"

using std::cout;
using std::endl;

Resonance::Resonance(int index, double E_cm, double dE_cm, double wg, double dwg, double Jr,
		     double G[3], double dG[3], int L[3], double Exf, int Int){
  // 'this' is a special pointer to the "current instance"
  this->index = index;
  this->E_cm = E_cm;
  this->dE_cm = dE_cm;
  this->wg = wg;
  this->dwg = dwg;
  this->Jr = Jr;
  for(int i=0; i<3; i++){
    this->G[i] = G[i];
    this->dG[i] = dG[i];
    this->L[i] = L[i];
  }
  this->Exf = Exf;
  this->bInt_flag = Int;

  //cout << "made a resonance" << endl;
}

Resonance::~Resonance(){}

void Resonance::print(){

  //  cout << "--------------------------------------------------" << "\n";
  //cout << "     This is resonance: " << index << "\n";
  cout << " Resonace " << std::setw(3) << index <<"    E_cm = " << E_cm  << " +/- " << dE_cm << "\n";
  cout << "                 wg   = " << wg << " +/- " << dwg << "\n";
  cout << "                 Jr   = " << Jr << "\n";
  cout << "                 G1   = " << G[0] << " +/- " << dG[0] << " (L = " << L[0] << ")\n";
  cout << "                 G2   = " << G[1] << " +/- " << dG[1] << " (L = " << L[1] << ")\n";
  cout << "                 G3   = " << G[2] << " +/- " << dG[2] << " (L = " << L[2] << ")\n";
  cout << "                 Exf  = " << Exf << "\n";
  cout << "           Integrated = " << bInt_flag << "\n";
  //  cout << "--------------------------------------------------" << "\n";
  cout << "\n";

}
