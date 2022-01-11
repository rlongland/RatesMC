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

#include "AMEReader.h"
#include "Utilities.h"


AMEReader::AMEReader(std::string filename) {

	// Open the file
	amefile.open(filename);

	if(!amefile.is_open()){
		std::cout << "ERROR: AME Mass file does not exist\n";
		exit(EXIT_FAILURE);			
	}
	
}

AMEReader::~AMEReader() {}

double AMEReader::readMass(std::string nucname) {

	std::cout << "The name is '" << nucname << "'\n";

	// Go through the file and look for where the name is the same as nucname
	std::string dummy;
	char name[2];
	char massfine[13];
	int massname;
	int masscourse;
	//	float massfine;
	int count=0;

	double mass=0.0;
	
	// Skip first lines
	for(int i=0; i<36; i++)
		std::getline(amefile,dummy);
	while(std::getline(amefile, dummy)){
		count++;
		sscanf(dummy.c_str(),"%*15c%3d %2s%*84c%3d %13c",&massname,name,&masscourse,massfine);
		std::stringstream readname;
		readname << massname << name;
		//		std::cout << "'" << readname.str() << "' " << masscourse << " " << massfine << std::endl;
		if(nucname.compare(readname.str()) == 0) break;
	}

	//	std::cout << std::scientific << std::setw(15) << std::setprecision(10) << masscourse << " " << massfine
	//						<< " " << masscourse + strtod(massfine,NULL)/1.0e6 << "\n";
	return masscourse + strtod(massfine,NULL)/1.0e6;
}
