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

void Reaction::getName(){
  std::cout << "The reaction name is: " << Name << std::endl;
}

void Reaction::printReaction(){

  cout << "--------------------------------------------------" << endl;
  cout << "     This is reaction: " << Name << endl;
  cout << "   Z0        Z1       Z2" << endl;
  printf( "   %2d        %2d       %2d\n",Z0,Z1,Z2);
  cout << "   M0        M1       M2" << endl;
  printf( "%5.2f     %5.2f    %5.2f\n",M0,M1,M2);
  cout << "   S_entrance = " << Q << endl;
  cout << "   S_exit     = " << Qexit << endl;
  cout << "   The gamma ray channel is channel " << Gamma_index << endl;
  cout << "--------------------------------------------------" << endl;

}
