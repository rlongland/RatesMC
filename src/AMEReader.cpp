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
/* ======================================================================
   AMEReader.cpp
   Description: Reads masses from the AME Mass table
   ======================================================================
*/
#include "stdio.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "AMEReader.h"
#include "Utilities.h"


AMEReader::AMEReader(std::string massfilename, std::string nubasefilename) {

	// Open the mass file
	amefile.open(massfilename);

	if(!amefile.is_open()){
		std::cout << "ERROR: AME Mass file does not exist\n";
		exit(EXIT_FAILURE);			
	}

	// Open the nubase file
	nubasefile.open(nubasefilename);

	if(!nubasefile.is_open()){
		std::cout << "ERROR: NuBase file does not exist\n";
		exit(EXIT_FAILURE);			
	}

}

AMEReader::~AMEReader() {}

double AMEReader::readMass(std::string nucname) {

	std::cout << "Reading AME mass for '" << nucname << "' = ";

	// Go through the file and look for where the name is the same as nucname

	std::string dummy;   // The AME line
	bool found = false;  // Was the nuclide found?
	
	double mass=0.0;     // The final mass

	int count = 0;
	// Seek back to the beginning
	amefile.seekg(0, amefile.beg);
	// Skip first 36 lines
	for(int i=0; i<36; i++)
		std::getline(amefile,dummy);
	// Now read line-by-line
	while(std::getline(amefile, dummy)){
		count++;

		std::string name;    // Element name
		int massname;        // Nuclide mass number (e.g. = 23 for 23Na)
		int masscourse;      // The integer part of the mass
		double massfine;     // The decimal part (in 10^-6 amu)

		// Use a stringstream to read the input
		std::stringstream ss(dummy);
		ss.ignore(16);// >> std::setw(16) >> dummy2;
		ss >> std::setw(3) >> massname;
		ss >> std::setw(3) >> name;
		ss.ignore(84);// >> std::setw(84) >> dummy2;
		ss >> std::setw(3) >> masscourse;
		ss >> std::setw(13) >> massfine;

		// Paste the mass number and element together. This is the nuclide's name
		std::stringstream readname;
		readname << massname << name;
		// Is this the nucleus we're looking for?
		if(nucname.compare(0,nucname.size(),readname.str()) == 0) {
			mass = masscourse + massfine/1.0e6;
			found = true;
			break;
		}
//		if(count > 300)break;
	}

	// Quit with an error if the nuclide was not found!
	if(found){
		std::cout << std::setprecision(8) <<  mass << " (atomic) ";
		return mass;
	} else {
		std::cout << "ERROR: Could not find AME mass for " << nucname << std::endl;
		exit(EXIT_FAILURE);
	}
}

//----------------------------------------------------------------------
double AMEReader::readMassFromAandZ(int A, int Z) {

	// Go through the file and look for where the name is the same as nucname
	std::cout << "Reading AME mass for A=" << A << " and Z=" << Z << " = ";
	
	std::string dummy;   // The AME line
	bool found = false;  // Was the nuclide found?
	
	double mass=0.0;     // The final mass

	int count = 0;
	// Seek back to the beginning
	amefile.seekg(0, amefile.beg);
	// Skip first 36 lines
	for(int i=0; i<36; i++)
		std::getline(amefile,dummy);
	// Now read line-by-line
	while(std::getline(amefile, dummy)){
		count++;

		std::string name;    // Element name
		int massname;        // Nuclide mass number (e.g. = 23 for 23Na)
		int chargename;      // Atomic number Z (e.g. = 11 for 23Na)
		int masscourse;      // The integer part of the mass
		double massfine;     // The decimal part (in 10^-6 amu)

		// Use a stringstream to read the input
		std::stringstream ss(dummy);
		ss.ignore(11);// >> std::setw(16) >> dummy2;
		ss >> std::setw(3) >> chargename;
		ss.ignore(2);
		ss >> std::setw(3) >> massname;
		ss >> std::setw(3) >> name;
		ss.ignore(84);// >> std::setw(84) >> dummy2;
		ss >> std::setw(3) >> masscourse;
		ss >> std::setw(13) >> massfine;

		// Is this the nucleus we're looking for?
		if(chargename == Z && massname == A) {
			mass = masscourse + massfine/1.0e6;
			found = true;
			break;
		}
//		if(count > 300)break;
	}

	// Quit with an error if the nuclide was not found!
	if(found){
		std::cout << std::setprecision(8) << mass << " (atomic) ";
		return mass;
	} else {
		std::cout << "ERROR: Could not find AME mass!" << std::endl;
		exit(EXIT_FAILURE);
	}
}


int AMEReader::readCharge(std::string nucname) {

	std::cout << "Reading AME charge for '" << nucname << "' = ";

	// Go through the file and look for where the name is the same as nucname

	std::string dummy;   // The AME line
	bool found = false;  // Was the nuclide found?
	
	int Z=0;     // The final mass

	int count = 0;
	// Seek back to the beginning
	amefile.seekg(0, amefile.beg);
	// Skip first 36 lines
	for(int i=0; i<36; i++)
		std::getline(amefile,dummy);
	// Now read line-by-line
	while(std::getline(amefile, dummy)){
		count++;

		std::string name;    // Element name
		int massname;        // Nuclide mass number (e.g. = 23 for 23Na)
		double charge;     // The atomic number (charge)

		// Use a stringstream to read the input
		std::stringstream ss(dummy);
		ss.ignore(11);// >> std::setw(16) >> dummy2;
		ss >> std::setw(3) >> charge;
		ss.ignore(2);
		ss >> std::setw(3) >> massname;
		ss >> std::setw(3) >> name;
		
		// Paste the mass number and element together. This is the nuclide's name
		std::stringstream readname;
		readname << massname << name;
		// Is this the nucleus we're looking for?
		if(nucname.compare(0,nucname.size(),readname.str()) == 0) {
			Z = charge;
			found = true;
			break;
		}
//		if(count > 300)break;
	}

	// Quit with an error if the nuclide was not found!
	if(found){
		std::cout << Z << "\n";
		return Z;
	} else {
		std::cout << "ERROR: Could not find AME charge for " << nucname << std::endl;
		exit(EXIT_FAILURE);
	}
}

//----------------------------------------------------------------------
// Read the spin 
double AMEReader::readSpin(std::string nucname) {

	std::cout << "Reading NuBase for '" << nucname << "': J = ";

	// Go through the file and look for where the name is the same as nucname

	std::string dummy;   // The AME line
	bool found = false;  // Was the nuclide found?
	
	double spin=0.0;     // The final mass

	int count = 0;
	// Seek back to the beginning
	nubasefile.seekg(0, nubasefile.beg);
	// Skip first 25 lines
	for(int i=0; i<25; i++)
		std::getline(nubasefile,dummy);
	// Now read line-by-line
	while(std::getline(nubasefile, dummy)){
		count++;

		std::string name;    // Element name
		std::string spinstring;
		double nucspin;     // The decimal part (in 10^-6 amu)

		// Use a stringstream to read the input
		std::stringstream ss(dummy);
		ss.ignore(11);// >> std::setw(16) >> dummy2;
		ss >> std::setw(5) >> name;
		ss.ignore(72);// >> std::setw(84) >> dummy2;
		ss >> std::setw(6) >> spinstring;

		// Get rid of any parentheses
		replace (spinstring.begin(), spinstring.end(), '(' , ' ');
		replace (spinstring.begin(), spinstring.end(), ')' , ' ');
		

		std::stringstream sj(spinstring);
		double numerator, denominator=1;
		char c;
		sj >> numerator;
		sj >> c;
		if(c == '/'){
			sj >> denominator;
		}
		nucspin = numerator/denominator;
		//std::cout << "'" << nucname << "' =? '" << name << "' J = '" << spinstring << "'\n";
		//std::cout << "Num = " << numerator << " den = " << denominator << "\n";
		// Is this the nucleus we're looking for?
		if(nucname.compare(0,nucname.size(),name) == 0) {
			spin = nucspin;
			found = true;
			break;
		}
		//if(count > 10)break;
	}

	// Quit with an error if the nuclide was not found!
	if(found){
		std::cout << spin << "\n";
		return spin;
	} else {
		std::cout << "ERROR: Could not find Nubase for " << nucname << std::endl;
		exit(EXIT_FAILURE);
	}
}
