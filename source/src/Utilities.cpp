#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <vector>
#include <sstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_coulomb.h>

#include "Utilities.h"

double EPS=1.0e-5;
double EMin;
gsl_rng * r;

std::ofstream logfile;
std::ofstream testfile;
std::ofstream ptfile;
int NSamples;
int NTemps;
bool ErrorFlag;
std::vector<double> Temp;

// counters
int PenZeroCount=0, IntegratedCount=0, SubSampledPosCount=0, SampledNegCount=0,
  NANCount=0, BelowIntLimit=0, IntfNANCount=0;

// Read an integer from a single line
int readInt(std::ifstream &infile){

  
  int input;
  infile >> input;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  return input;

}
// Read a double from a single line
double readDouble(std::ifstream &infile){

  double input;
  infile >> input;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  return input;

}
void skipLines(std::ifstream &infile, int nlines){

  std::string dummy;
  for(int i=0; i<nlines; i++){
    infile >> dummy;
    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
  
}

// Count the number of entries on a line
int countEntries(std::ifstream &infile){

  std::string line;
  std::vector<std::string> entries;

  // First save position
  int place=infile.tellg();

  std::getline(infile, line);
  std::stringstream ss(line);
  std::string entry;
  while( ss >> entry ){
    entries.push_back(entry);
    //std::cout << entry << " ";
  }
  //  std::cout << "\n";

  // Return to saved position in file
  infile.seekg(place);
  
  return entries.size();
  
}

void readNonResonant(std::ifstream &infile, Reaction &R, int part){

  int ne = countEntries(infile);
  //std::cout << "There are " << ne << " entries\n";
  if(ne != 5){
    std::cout << "  ERROR: There should be 5 numbers in the non-resonant input lines\n";
    exit(EXIT_FAILURE);
  }

  double s, sp, spp, ds, cutoffe;
  infile >> s >> sp >> spp >> ds >> cutoffe;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R.setNonResonant(s,sp,spp,ds,cutoffe,part);

}

/* 
   This is where Victor need to work on reading all of the "standard" resonances
   1) Start by just reading 5 of them
   2) Then add the ability to add up until the end of the block (check for ****?)
   3) Include error checking. Do the numbers make sense?
*/

void readResonanceBlock(std::ifstream &infile, Reaction &R, bool isUpperLimit){

  double E_cm, dE_cm, wg=0.0, dwg=0.0, Jr,
    G1, dG1, PT1=0.0, DPT1=0.0, G2, dG2, PT2=0.0, DPT2=0.0,
    G3, dG3, PT3=0.0, DPT3=0.0, Exf;
  int i,L1, L2, L3, isBroad;

  i=0;
  // Read the number of entries on the first resonance line. This
  // gives us another check on whether it's an upper limit or normal
  // resonance
  int nEnt = countEntries(infile);
  //  std::cout << "Reading resonances:\n";
  //std::cout << "There are " << countEntries(infile) << " entries in this section\n";
  if(!isUpperLimit){
    isUpperLimit = false;
    //std::cout << "Normal resonances\n";
    if(nEnt != 16){
      std::cout << "ERROR! The number of columns in the resonance section is wrong\n";
      std::cout << "       Expect N=16; Got N=" << nEnt << "\n";
    }
  } else {
    isUpperLimit = true;
    //std::cout << "Upper limit resonances\n";
    if(nEnt == 17){
      std::cout << "WARNING: It looks like you're using an old version of RatesMC input without DPT\n";
    } else if (nEnt == 20) {
    } else {
      std::cout << "ERROR! The number of columns in the upper limit resonance section is wrong\n";
      std::cout << "       Expect N=17 or 20; Got N=" << nEnt << "\n";
    }
  }
     
  while(true){

    // First try to read resonance energy to see if it's a real resonance input
    std::string data;
    infile >> data;
    /*
    std::stringstream ss(data);
    std::string entry;
    ss >> entry; 
    if(!std::isdigit( entry[0] ) && !(entry[0] == '-') ) {
      //      std::cout << "Found end of resonances!\n" << data << "\n";
      break;
    }
    */
    // Look for '***', which signifies the end of resonance input
    std::size_t found;
    found = data.find("***");
    if (found!=std::string::npos){
      logfile << "Found end of ";
      if(isUpperLimit)
	logfile << "upper limit ";
      logfile << "resonances!\n" << data << std::endl;
      break;
    }

    // Now look for an exclamation point, which indicates a commented-out resonance
    found = data.find("!");
    if (found!=std::string::npos){
      if(isUpperLimit)
	logfile << "Upper limit ";
      logfile << "Resonance is commented-out! " << data << std::endl;;
      continue;
    }
    
    // If this looks like a resonance, read a single line
    E_cm = std::stod(data);
    if(!isUpperLimit){
      infile >> dE_cm >> wg >> dwg >> Jr 
	     >> G1 >> dG1 >> L1 >> G2 >> dG2 >> L2 >> G3 >> dG3 >> L3
	     >> Exf >> isBroad;
    } else {
      // Upper limit resonances.
      // Old style with no DPT
      if(nEnt == 17){
	infile >> dE_cm >> Jr 
	       >> G1 >> dG1 >> L1 >> PT1 >> G2 >> dG2 >> L2 >> PT2
	       >> G3 >> dG3 >> L3 >> PT3
	       >> Exf >> isBroad;
	
      }
      // New style with DPT
      else if(nEnt == 20){
	infile >> dE_cm >> Jr 
	       >> G1 >> dG1 >> L1 >> PT1 >> DPT1 >> G2 >> dG2 >> L2 >> PT2 >> DPT2
	       >> G3 >> dG3 >> L3 >> PT3 >> DPT3
	       >> Exf >> isBroad;
      }
    }

    // Convert to correct units
    E_cm  *= 1.0e-3;   // to be in MeV
    dE_cm *= 1.0e-3;   // to be in MeV
    wg    *= 1.0e-6;   // to be in eV
    dwg   *= 1.0e-6;   // eV
    G1    *= 1.0e-6;   // eV
    dG1   *= 1.0e-6;   // eV
    G2    *= 1.0e-6;   // eV
    dG2   *= 1.0e-6;   // eV
    G3    *= 1.0e-6;   // eV
    dG3   *= 1.0e-6;   // eV

    infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

    // Add this resonance to the list of resonances stored in the
    // reaction
    R.addResonance(i++, E_cm, dE_cm, wg, dwg, Jr,
		   G1, dG1, L1, PT1, DPT1,
		   G2, dG2, L2, PT2, DPT2,
		   G3, dG3, L3, PT3, DPT3,
		   Exf, isBroad, isUpperLimit);
    
  }

  
  
}

int ReadInputFile(std::string inputfilename, Reaction *R){

  std::cout << "The input file name is: " << inputfilename << std::endl;

  int ne;  // for counting entries
  
  // Start by opening input file (reactions.dat)
  std::ifstream infile;
  infile.open(inputfilename.c_str(),std::ifstream::in);

  // test to see if file is there
  if(!infile.is_open()){
    std::cout << "ERROR! File: RatesMC.in is not found!" << std::endl;
    return(1);
  }

  std::string name,dummy;
  infile >> name;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setName(name);
  //R -> getName();

  // Blank line
  infile >> dummy;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  // Masses, charges, spins, Q-values at the top
  int z0,z1,z2;
  double R0;
  double m0,m1,m2;
  double j0,j1,j2;
  double Qin,Qout;
  int gindex;
  
  infile >> z0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> z1;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> z2;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setCharges(z0,z1,z2);

  infile >> m0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> m1;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> m2;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setMasses(m0,m1,m2);

  infile >> j0;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> j1;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  infile >> j2;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  R -> setSpins(j0,j1,j2);

  // Entrance and exit particle separation energies
  Qin = readDouble(infile);
  Qout = readDouble(infile);
  R -> setSeparationEnergies(Qin, Qout);

  // The R0 radius
  R0 = readDouble(infile);
  R -> setR0(R0);

  // Index of the gamma-ray channel
  gindex = readInt(infile);
  R -> setGammaIndex(gindex);

  // Ignore a line
  skipLines(infile, 1);

  // Computation control block
  EMin = readDouble(infile);
  
  NSamples = readInt(infile);
  NTemps = readInt(infile);
  int itmp;
  infile >> itmp;
  bool bPartialWidthCorrelations = (bool)itmp;
  infile >> itmp;
  bool bEnergyCorrelations = (bool)itmp;
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

  // Skip 3 lines
  skipLines(infile, 3);

  // Non-resonant line 1
  readNonResonant(infile, *R, 0);
  
  // Non-resonant line 2
  readNonResonant(infile, *R, 1);

  // Skip 5 lines
  skipLines(infile, 5);

  // Read all of the resonances 
  readResonanceBlock(infile, *R, false);

  // Skip 4 lines
  skipLines(infile, 4);
  // Read the upper limit resonances
  readResonanceBlock(infile, *R, true);
  
  return 0;
}

//----------------------------------------------------------------------
// Define the temperature array
void defineTemperatures(){

  std::vector<double> defaultT{0.01,0.011,0.012,0.013,0.014,0.015,
			 0.016,0.018,0.020,0.025,0.03,0.04,
			 0.05,0.06,0.07,0.08,0.09,0.1,0.11,
			 0.12,0.13,0.14,0.15,0.16,0.18,0.20,
			 0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,
			 0.8,0.9,1.0,1.25,1.5,1.75,2,2.5,3,
			 3.5,4,5,6,7,8,9,10};
  
  Temp = defaultT;

  logfile << "--------------------------------------------------\n";
  logfile << "There are " << Temp.size() << " temperatures:\n";
  for(int iT=0; iT<Temp.size(); iT++){
    logfile << std::setw(5) << Temp[iT] << "  ";
    if((iT+1)%7 == 0 && iT>0)logfile << "\n";
  }
  logfile << "\n";
  logfile << "--------------------------------------------------";
  logfile << std::endl;

}

//----------------------------------------------------------------------
void logNormalize(double mean, double sd, double& mu, double& sigma){

  // Varience is standard deviation squared
  double var = gsl_pow_2(sd);

  if(mean == 0.0){
    mu = 0.0;
    sigma = 0.0;
    return;
  }

  mu = gsl_sf_log(mean) - 0.5*gsl_sf_log(1+(var/gsl_pow_2(mean)));
  sigma = sqrt(gsl_sf_log(1+(var/gsl_pow_2(mean))));

  // If mean and var are very small, some computers will return nan for these ratios
  // so check for it
  if(gsl_isnan(mu) || gsl_isnan(sigma)){
    mu = 0.0;
    sigma = 0.0;
  }

  return;


}


//----------------------------------------------------------------------
// Penetration factor calculation
double PenFactor(double E, double L, double Mass0, double Mass1,
		 int Charge0, int Charge1, double R){

  // Turn off the GSL error handler, which aborts the program
  // if G or F go out of range.
  gsl_set_error_handler_off();

  gsl_sf_result F,Fp,G,Gp;
  double exp_F,exp_G,P;
  double mue = Mass0*Mass1/(Mass0+Mass1);

  //  double R = R0*(pow(Mass0,(1./3.))+pow(Mass1,(1./3.)));
  double rho = 0.218735097*R*sqrt(mue*E);

  double eta = 0.15748927*(double)Charge0*(double)Charge1*sqrt(mue/E);

  int status = gsl_sf_coulomb_wave_FG_e (eta, rho, L, 0, &F, &Fp, &G,
					 &Gp, &exp_F, &exp_G);

  // Check to see if this failed. If it's out of range, set P=0,
  // if not, print an error and exit
  if(status){
    if(status == GSL_EOVRFLW){
      ErrorFlag = true;
      PenZeroCount++;
      return 0.0;
    } else {
      std::cout << "\nERROR: Something went wrong in coulomb wavefunction!" <<
	"\n\tGSL Error: " << gsl_strerror (status)<< std::endl;
      std::cout << "The Energy was " << E*1e3 << " keV." << std::endl;
      abort();
    }
  }
  gsl_set_error_handler (NULL);


  // Just in case there is an overflow, multiply by exponential
  // (See GSL documentation for more info)
  double F_l = F.val*exp(exp_F);
  double G_l = G.val*exp(exp_G);

  P = rho/(gsl_pow_2(F_l) + gsl_pow_2(G_l));

  return P;
}

//----------------------------------------------------------------------
// Check for zero
bool isZero(double x){
  const double epsilon = 1e-5;
  return std::abs(x-0.) <= epsilon*std::abs(x);
}

//----------------------------------------------------------------------
// Setup the random sampler
void setupRandom(){
  
   // set-up the GSL random sampler
  const gsl_rng_type * T;
  /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  // Get a seed from either /dev/random, or the clock
  unsigned long int seed = random_seed();
  // Set the seed in the RNG
  gsl_rng_set(r,seed);

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Function to find the random seed
// Courtesy of Robert G. Brown

#include <sys/time.h>

unsigned long int random_seed()
{

  unsigned int seed;
  int iret=0;
  struct timeval tv;
  FILE *devrandom;

  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
    //   if(verbose == D_SEED) printf("Got seed %u from gettimeofday()\n",seed);
  } else {
    iret = fread(&seed,sizeof(seed),1,devrandom);
    //   if(verbose == D_SEED) printf("Got seed %u from /dev/random\n",seed);
    fclose(devrandom);
  }
  if(iret < 0)std::cout << "WARNING! There could be a problem with the random numbers!\n";
  return(seed);

}
