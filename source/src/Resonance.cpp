/* ======================================================================
   Resonance.cpp
   Author: R. Longland
   Date: 2021-01-27
   
   Description: Contains all of the resonance specific stuff
   ======================================================================
*/
#include <iostream>
#include "stdio.h"
#include "Resonance.h"

using std::cout;
using std::endl;

Resonance::Resonance(int index, double E_cm, double dE_cm, double wg, double dwg){
  // 'this' is a special pointer to the "current instance"
  this->index = index;
  this->E_cm = E_cm;
  this->dE_cm = dE_cm;
  this->wg = wg;
  this->dwg = dwg;
  cout << "made a resonance" << endl;
}

Resonance::~Resonance(){}

void Resonance::print(){

  cout << "--------------------------------------------------" << "\n";
  cout << "     This is resonance: " << index << "\n";
  cout << "                 E_cm = " << E_cm  << "\n";
  cout << "--------------------------------------------------" << "\n";
  cout << "\n";

}
