/* ======================================================================
   Reaction.cpp
   Author: R. Longland
   Date: 2018-11-13
   
   Description: Contains all of the reaction specific stuff
   ======================================================================
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "stdio.h"
#include "Resonance.h"
#include "Reaction.h"
#include "Utilities.h"

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

void Reaction::addResonance(int i, double E_cm, double dE_cm, double wg, double dwg, double Jr,
			    double G1, double dG1, int L1, double G2, double dG2, int L2,
			    double G3, double dG3, int L3, double Exf, int Int){

  double G[3] = {G1, G2, G3};
  double dG[3] = {dG1, dG2, dG3};
  int L[3] = {L1, L2, L3};
  
  // Make a resonance
  Resonance Res(i, E_cm, dE_cm, wg, dwg, Jr,
		G, dG, L, Exf, Int);
  //Res.print();
  // Add that resonance to the list of resonances
  Resonances.push_back(Res);
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

  cout << "--------------------------------------------------" << "\n";
  cout << " Resonances:\n";
  // Loop through all regular resonances
  std::vector<Resonance>::iterator res;
  for(res = Resonances.begin(); res < Resonances.end(); res++){
    res->print();
  }

}

void Reaction::writeReaction(){

  logfile << "--------------------------------------------------" << "\n";
  logfile << "     This is reaction: " << Name << "\n";
  logfile << "   Z0        Z1       Z2" << "\n";
  logfile << std::setw(5) << Z0 << "     " << std::setw(5) << Z1 << "    " << std::setw(5) << Z2 << "\n";
  logfile << "   M0        M1       M2" << "\n";
  logfile << M0 << "   " << M1 << "   " << M2 << "\n";
  logfile << "   S_entrance = " << Q << "\n";
  logfile << "   S_exit     = " << Qexit << "\n";
  logfile << "   The gamma ray channel is channel " << Gamma_index << "\n";
  logfile << "--------------------------------------------------" << "\n";

  logfile << " Direct Capture part     \n";
  logfile << "          S0       S'       S''    dS   CutoffE\n";
  logfile << " Part 1: " << S[0] << "  " << Sp[0] << "  " << Spp[0] << "  " << dS[0] << "  " << CutoffE[0] << "\n";
  logfile << " Part 2: " << S[1] << "  " << Sp[1] << "  " << Spp[1] << "  " << dS[1] << "  " << CutoffE[1] << "\n";
  logfile << "\n";

  logfile << "--------------------------------------------------" << "\n";
  logfile << " Resonances:\n";
  // Loop through all regular resonances
  std::vector<Resonance>::iterator res;
  for(res = Resonances.begin(); res < Resonances.end(); res++){
    res->write();
  }

}
