/* RatesMC: Monte Carlo reaction rate code
 *
 * Copyright (C) 2022  Richard Longland
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * Author: R. Longland
 *  Email: rllongla@ncsu.edu
 *
 */

#ifndef _AMEReader_h_
#define _AMEReader_h_

#include <iostream>

class AMEReader {

 public:

  // Constructor
  AMEReader(std::string massfilename, std::string nubasefilename);
  // Destructor
  ~AMEReader();

  double readMass(std::string nucname);
  double readMassFromAandZ(int A, int Z);
	int readCharge(std::string nucname);
  double readSpin(std::string nucname);

	void toggleAMEcharge(int i){AMEcharge[i]=true;}
	void toggleAMEmass(int i){AMEmass[i]=true;}

	bool getAMEmass(int i){return AMEmass[i];}
	
 private:

  std::ifstream amefile;
	std::ifstream nubasefile;

	// flags
	bool AMEcharge[3]={false};
	bool AMEmass[3]={false};
	
};

#endif
