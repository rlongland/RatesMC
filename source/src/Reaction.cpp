/* ======================================================================
   Reaction.cpp
   Author: R. Longland
   Date: 2018-11-13
   
   Description: Contains all of the reaction specific stuff
   ======================================================================
*/
#include <iostream>
#include "stdio.h"
#include "Reaction.h"

using std::cout;
using std::endl;

Reaction::Reaction(){}
Reaction::~Reaction(){}

void Reaction::setNonResonant(double s, double sp, double spp, double ds, double cutoffe, int part){
  S[part] = s;
  Sp[part] = sp;
  Spp[part] = spp;
  dS[part] = ds;
  CutoffE[part] = cutoffe;
}


void Reaction::getName(){
  cout << "The reaction name is: " << Name << "\n";
}

void Reaction::printReaction(){

  cout << "--------------------------------------------------" << "\n";
  cout << "     This is reaction: " << Name << "\n";
  cout << "   Z0        Z1       Z2" << "\n";
  printf( "   %2d        %2d       %2d\n",Z0,Z1,Z2);
  cout << "   M0        M1       M2" << "\n";
  printf( "%5.2f     %5.2f    %5.2f\n",M0,M1,M2);
  cout << "   S_entrance = " << Q << "\n";
  cout << "   S_exit     = " << Qexit << "\n";
  cout << "   The gamma ray channel is channel " << Gamma_index << "\n";
  cout << "--------------------------------------------------" << "\n";

  cout << " Direct Capture part     \n";
  cout << "          S0       S'       S''    dS   CutoffE\n";
  cout << " Part 1: " << S[0] << "  " << Sp[0] << "  " << Spp[0] << "  " << dS[0] << "  " << CutoffE[0] << "\n";
  cout << " Part 2: " << S[1] << "  " << Sp[1] << "  " << Spp[1] << "  " << dS[1] << "  " << CutoffE[1] << "\n";
  cout << "\n";

}
