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
   Interference.cpp
   Description: Container to hold an interfering pair of resonances
   ======================================================================
*/
#include "stdio.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <limits>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

#include "Interference.h"
#include "Resonance.h"
#include "Utilities.h"

// Ugly-ass hack
Interference *InterferencePtr;
double InterferenceIntegrandWrapper(double x, void *params) {
  return InterferencePtr->Integrand(x, params);
}
// end ugly-ass hack
// another ugly-ass hack
int rhsInterferenceWrapper(double x, const double y[], double dydt[], void *params_ptr){
	InterferencePtr->rhs(x, y, dydt, params_ptr);
	return GSL_SUCCESS;
}
// end

Interference::Interference(Reaction &R, int i, int isign)
	: Reac(R) {

  index = i;
	IntfSign = isign;

}

Interference::~Interference() {}

void Interference::addResonance(Resonance* R, int index) {
	std::cout << "Adding resonance to part " << index << "\n"; 
	Res[index] = R;
	
}

void Interference::makeSamples(std::vector<std::vector<double> > Ref_sample,
															 double smallestdE,
															 double smallestdwg, double smallestdG[3]){

	// Make sure the Rate vector is the right size
  Rate_sample.resize(NSamples);

	
	Res[0]->makeSamples(Ref_sample, smallestdE, smallestdwg, smallestdG);
	Res[1]->makeSamples(Ref_sample, smallestdE, smallestdwg, smallestdG);

	M0 = Res[0]->getM0();
	M1 = Res[0]->getM1();
	M2 = Res[0]->getM2();
	Z0 = Res[0]->getZ0();
	Z1 = Res[0]->getZ1();
	Z2 = Res[0]->getZ2();
	J0 = Reac.J0;
	J1 = Reac.J1;
	
  R = Reac.R0 * (pow(M0, (1. / 3.)) + pow(M1, (1. / 3.)));
  mue = M0 * M1 / (M0 + M1);

	for(int iRes = 0; iRes<2; iRes++){
		
		E_cm[iRes] = Res[iRes]->getE_cm();
		Jr[iRes] = Res[iRes]->getJr();
		for(int i=0; i<3; i++){
			G[iRes][i] = Res[iRes]->getG(i);
			L[iRes][i] = Res[iRes]->getL(i);
		}
		Exf[iRes] = Res[iRes]->getExf();
		NChannels[iRes] = Res[iRes]->getNChannels();
		//std::cout << "MakeSamples: Exf[" << iRes << "] = " << Res[iRes]->getExf() << std::endl;
	}
	
	
	// get the interference sign
	switch(IntfSign){
	case -1:
		sign_sample.resize(NSamples);
		std::fill(sign_sample.begin(),sign_sample.end(),-1);
		break;
	case 1:
		sign_sample.resize(NSamples);
		std::fill(sign_sample.begin(),sign_sample.end(),1);
		break;
	case 0:
		for(int i=0; i<NSamples; i++)
			sign_sample.push_back(2*gsl_rng_uniform_int (r,2)-1);
		break;
	}

	
}

//----------------------------------------------------------------------
// Function to numerically integrate broad resonances
double Interference::calcBroad(double T) {

  double classicalRate = 0.0; //, ARate;


  // Calculate the rate samples
  for (int s = 0; s < NSamples; s++) {

		if (s % 10 == 0) {
			// \r goes back to the beginning of the line.
			std::cout << "\r" << 100 * s / NSamples << "% Complete for Interference "
								<< index + 1 << std::flush;
		}

		// Calculate the single integrated rate sample
		// For each resonance we need to grab the samples
		double E_sample[2];
		double G_sample[2][3];
		double erFrac[2][3];
		
		for(int iRes=0; iRes<2; iRes++){
			E_sample[iRes] = Res[iRes]->getESample(s);
			for(int i=0; i<3; i++){
				G_sample[iRes][i] = Res[iRes]->getGSample(i,s);
				erFrac[iRes][i] = Res[iRes]->geterFrac(i,s);
			}			
		}
    Rate_sample[s] = NumericalRate(T, E_sample, G_sample, erFrac, sign_sample[s], false);
		//std::cout << "Returned from NumericalRate with " << Rate_sample[s] << "\n";
  }

  // And the central value, which is the classical rate
	// Write the integrand to a file
	double erFrac[2][3];
	for(int iRes=0; iRes<2;iRes++)
		for(int i=0; i<3; i++)
			erFrac[iRes][i] = 1.0;
	switch(IntfSign){
	case(0):
		classicalRate = NumericalRate(T, E_cm, G, erFrac, 1, true);
		classicalRate = NumericalRate(T, E_cm, G, erFrac, -1, true);
		break;
	case(1):
		classicalRate = NumericalRate(T, E_cm, G, erFrac, 1, true);
		break;
	case(-1):
		classicalRate = NumericalRate(T, E_cm, G, erFrac, -1, true);
		break;
	}
	integrandfile << std::endl;

	//	std::cout << "ClassicalRate = " << classicalRate << "\n";
	
  return classicalRate;
}

//----------------------------------------------------------------------
// Multiple the rate samples by a constant
void Interference::scaleByConstant(double RateFactor) {
  std::transform(Rate_sample.begin(), Rate_sample.end(), Rate_sample.begin(),
                 [RateFactor](double &c) { return c * RateFactor; });
}
//----------------------------------------------------------------------
// Print out the rate
void Interference::printRate() {

  // Convert all vector to an array
  double *Rate_Array = &Rate_sample[0];

  // If there are enough samples, calculate the mean, variance,
  // log-mean and log-variance and bin the rates
  // Turn off the gsl error handler in case we have a negative number here
  gsl_set_error_handler_off();
  double lograte[NSamples];
  gsl_sf_result LogResult;
  if (NSamples > 2) {
    for (int s = 0; s < NSamples; s++) {

      int status = gsl_sf_log_e(Rate_Array[s], &LogResult);
      // Check to make sure log didn't error
      if (status) {
        // ErrorFlag = true;
        // LogZeroCount++;
        lograte[s] = 0.0;
      } else {
        lograte[s] = LogResult.val;
      }
    }
    MeanRate = gsl_stats_mean(Rate_Array, 1, NSamples);
    MedianRate = gsl_stats_median(Rate_Array, 1, NSamples);
    RateMu = gsl_stats_mean(lograte, 1, NSamples);
    RateSigma = sqrt(gsl_stats_variance(lograte, 1, NSamples));
  }

  std::cout << "Mean = " << MeanRate << " Median = " << MedianRate
            << " Mu = " << RateMu << " Sigma = " << std::scientific << RateSigma
            << "\n";
}

//----------------------------------------------------------------------
// This function gets called for each sample
double Interference::NumericalRate(double T, double E_sample[2],
																	 double G_sample[2][3],
																	 double erFrac[2][3],
																	 double sign_sample,
																	 bool writeIntegrand) {

  double ARate = 0.0;
		
	double Pr[2], Pr_exit[2];

  InterferencePtr = this;

  // if reaction is endothermic, need to make minimum energy
  // enough to ensure no integration over negative energies
  double E_min = EMin;
	double E_max = 10.0;

	// For each resonance in the interfering pair, do any
	// energy-dependent scaling, sub-threshold resonance conversion,
	// etc.
	for(int iRes = 0; iRes<2; iRes++){
		//		std::cout << "iRes = " << iRes << " with E = " << E_cm[iRes] << "\n";
		//std::cout << "E_sample = " << E_sample[iRes] << ", L[iRes][] = " << L[iRes][0] << " " <<
		//	L[iRes][1] << " " << L[iRes][2] << "\n";
		// if particle is in spectator channel, integration should
		//  not be truncated
		if (Reac.Qexit > Reac.Q && Reac.getGamma_index() == 2)
			E_min += Reac.Qexit + Exf[iRes] - Reac.Q;

		// Scale the partial widths by the energy effect
		for(int i=0; i<3; i++){
			G_sample[iRes][i] *= erFrac[iRes][i];
		}
	
		// IF we input E>0, but the sample is <0, we need to treat it as a
		// subthreshold resonance. To do this, we need to convert
		// G into C2S*Theta_sp
		//std::cout << "here\n";
		if (E_cm[iRes] > 0.0 && E_sample[iRes] < 0.0) {
			ErrorFlag = true;
			SampledNegCount++;
			G_sample[iRes][0] = mue * gsl_pow_2(R) * G_sample[iRes][0] /
				(2.0 * 41.80161396 * PenFactor(E_cm[iRes], L[iRes][0], M0, M1, Z0, Z1, R));
			//std::cout << "Positive resonance went negative!\n";
			//std::cout << "E_cm = " << E_cm << "E_sample = " << E_sample[iRes] << "G[0] = " <<
			//	G[0]	<< " G_sample = " << G_sample[iRes][0] << std::endl;
		} else if (E_cm[iRes] < 0.0 && E_sample[iRes] > 0.0) {
			// or convert to real resonance if E>0.0
			ErrorFlag = true;
			SubSampledPosCount++;
			G_sample[iRes][0] = G_sample[iRes][0] * 2.0 * 41.80161396 *
				PenFactor(E_sample[iRes], L[iRes][0], M0, M1, Z0, Z1, R) /
				(mue * gsl_pow_2(R));
			//	std::cout << "Negative resonance went positive!\n";
			//		std::cout << "E_cm = " << E_cm << "E_sample = " << E_sample[iRes] << "G[0] = " << G[0]
			//							<< " G_sample = " << G_sample[iRes][0] << std::endl;
		}

		//  The penetration factor at the resonance energy (the "true" PF)
		//std::cout << "normal\n";
		if (E_sample[iRes] > 0.0) {
			//std::cout << M0 << " " << M1 << " " << M2 << " " << Z0 << " " << Z1 << " " << R << "\n";
			Pr[iRes] = PenFactor(E_sample[iRes], L[iRes][0], M0, M1, Z0, Z1, R);
			//std::cout << "Pr = " << Pr[iRes] << "\n";
			//if(isZero(Pr[iRes]))return 0.0;
		} else {
			Pr[iRes] = 0.0;
		}
		//		std::cout << "Pr = " << Pr[iRes] << "\n";

	
		// Calculate the exit particle energy, depends on if it is spectator
		// if(NChannels[j]==3){
		if (Reac.getGamma_index() == 2) {
			// if exit particle is observed decay, take final excitation into account
			Pr_exit[iRes] = PenFactor(E_sample[iRes] + Reac.Q - Reac.Qexit - Exf[iRes], L[iRes][1],
																M0 + M1 - M2, M2,
																Z0 + Z1 - Z2, Z2, R);
			// cout << "Exit energy = " << E+Reac.Q-Reac.Qexit-Exf << endl;
		} else if (Reac.getGamma_index() == 1 && NChannels[iRes] == 3) {
			// ignore spectator final excitation if it is spectator
			Pr_exit[iRes] = PenFactor(E_sample[iRes] + Reac.Q - Reac.Qexit, L[iRes][2],
																M0 + M1 - M2, M2,
																Z0 + Z1 - Z2, Z2, R);
			// cout << "Exit energy = " << E+Reac.Q-Reac.Qexit << endl;
		}
	}

	// Now we're ready to do the integration!
		
	//--------------------------------------------------
	// GSL Integration functions
  double result, error;
  size_t nevals;

  // Define integration limits
  //double x = E_min, x1 = E_max;
  //double pole = E;

  // The array that needs to be passed to the integration function
	//	std::cout << "About to start integration\n";
  double alpha[15];
  alpha[0] = (int)writeIntegrand;
  alpha[1] = Pr[0];
  alpha[2] = Pr[1];
  alpha[3] = Pr_exit[0];
  alpha[4] = Pr_exit[1];
  alpha[5] = E_sample[0];
  alpha[6] = E_sample[1];
  alpha[7] = G_sample[0][0];
  alpha[8] = G_sample[1][0];
  alpha[9] = G_sample[0][1];
  alpha[10] = G_sample[1][1];
  alpha[11] = G_sample[0][2];
  alpha[12] = G_sample[1][2];
  alpha[13] = sign_sample;
  alpha[14] = T;

  //double gammaT = G0 + G1 + G2;

  //  std::cout << "making integration workspace\n";
  // alpha[0] = 1;

  gsl_function F;
  // Can't use Integrand directly because GSL is shit
  F.function = &InterferenceIntegrandWrapper;
  //  F.function = &Integrand;
  F.params = &alpha;

  /*
		// Some integration methods that I had trouble with....
   int status = gsl_integration_qags (&F,      // Function to be integrated
                        E_min,       // start
                        E,      // end
                        0,       // absolute error
                        1e-3,    // relative error
                        1000,    // max number of steps (cannot exceed size of
  workspace w,       // workspace &result, // The result &error); ARate =
  result; status = gsl_integration_qags (&F,      // Function to be integrated
                        E,       // start
                        E_max,      // end
                        0,       // absolute error
                        1e-3,    // relative error
                        1000,    // max number of steps (cannot exceed size of
  workspace w,       // workspace &result, // The result &error); ARate +=
  result;
  */


  /*
  double pts[3];
  pts[0] = E_min;
  pts[1] = E;
  pts[2] = E_max;

  // If it's subthreshold
  if(E < 0.0){npts
    pts[1] = pts[0];
  }
  */

  //  std::cout << "Integration pts = " ;
  // for(int i=0; i<npts; i++) std::cout << pts[i] << " ";
  // std::cout << std::endl;

  // if(writeIntegrand)
  //  std::cout << E_min << " " << gammaT << " " << E << " " << pts[2]-pts[1] <<
  //  " " << E_max << "\n";

  //        std::cout << E << "\n";

	// Finally settles on cquad integration routine
  // Turn off the error handler
  gsl_error_handler_t *temp_handler;
  temp_handler = gsl_set_error_handler_off();

	//std::cout << "Setting up integration workspace\n";
  gsl_integration_cquad_workspace *w =
      gsl_integration_cquad_workspace_alloc(10000);
	//std::cout << "Integration!\n";
  int status = gsl_integration_cquad(&F,      // Function to be integrated
                                     E_min,   // Where known singularity is
                                     E_max,   // number of singularities
                                     1e-100,   // absolute error
                                     1e-6,    // relative error
                                     w,       // workspace
                                     &result, // The result
                                     &error, &nevals);
  gsl_integration_cquad_workspace_free(w);

	//std::cout << "Done integration: " << status << "\n";
	//status = -1;
  // If the integration errored, use the slower ODE method
  if (status != 0) {
		//std::cout << "Integration error = " << gsl_strerror(status) << "\n";
    //result = std::numeric_limits<double>::quiet_NaN();
	
		// OK so the fast integration failed. Go back to the old method
		const gsl_odeiv2_step_type * T
			= gsl_odeiv2_step_rkck;

		// Set up the ODE solver
		gsl_odeiv2_step * s
			= gsl_odeiv2_step_alloc (T, 1);
		gsl_odeiv2_control * c
			= gsl_odeiv2_control_y_new (1e-100, 1.0e-6);
		gsl_odeiv2_evolve * e
			= gsl_odeiv2_evolve_alloc (1);

		// Function, Jacobian, Number of Dimensions, Parameters
		gsl_odeiv2_system sys = {rhsInterferenceWrapper, NULL, 1, &alpha};

		// Integration limits
		double x = E_min, x1 = E_max;
		// stepsize
		double h=1.0e-10;
		double hmin = 1.0e-12;   // The minimum step size
		// Value of the integrand. Starts at zero
		double y[2] = {0.0, 0.0};
		// Flag to take small steps
		bool smallstep = false;

		while(x < x1){

      // 2007-12-27
      // If the step will bring it close to the Er, take small steps
      if(x < E_sample[0] && (x + h) > (E_sample[0]-(3.0*(G_sample[0][0]+
																												 G_sample[0][1]+
																												 G_sample[0][2]))))smallstep=true;
      if(x < E_sample[1] && (x + h) > (E_sample[1]-(3.0*(G_sample[1][0]+
																												 G_sample[1][1]+
																												 G_sample[1][2]))))smallstep=true;
			
      // Now make step small if it isn't already, can be fairly large because
      // it starts on the resonance wing
      if(smallstep && h>(G_sample[0][0]+G_sample[0][1]+G_sample[0][2])){
				h = (G_sample[0][0]+G_sample[0][1]+G_sample[0][2]);
      }
      if(smallstep && h>(G_sample[1][0]+G_sample[1][1]+G_sample[1][2])){
				h = (G_sample[1][0]+G_sample[1][1]+G_sample[1][2]);
      }

      // Make the step...
      // Note, the step size is determined by the error in the
      //  cross section, so spectroscopic factors may look jumpy
      int status = gsl_odeiv2_evolve_apply(e, c, s,
					  &sys,
					  &x, x1,
					  &h, y);

      // quit if there was an error
      if(status != GSL_SUCCESS)break;

      // don't let h get too small
      if(h<hmin)h=hmin;

      // reset small step flag
      smallstep=false;
      //     if(E < 0.05)testhist << x << "\t" << y[0] << endl;
      //cout << h << "\t" << y[0] << endl;
      //testhist << x << "\t" << y[0] << endl;

    }

		result = y[0];

	}
  gsl_set_error_handler(temp_handler);

	//	std::cout << "Result = " << result << "\n";
	// The integration result
  ARate = result;

  // If G0 or G1 were zero, sum is NAN, catch this!
  if (isnan(ARate)) {
    // ARate = 0.0;
    ErrorFlag = true;
    NANCount++;
  }
  if (isinf(ARate)) {
    // Arate = 0.0;
    ARate = std::numeric_limits<double>::quiet_NaN();
    ErrorFlag = true;
    InfCount++;
  }

  return ARate;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Function to be integrated
// extern "C" {
double Interference::Integrand(double x, void *params) {

  //  std::cout << "Integrand\n";

  double *par = (double *)params;
  int writeIntegrand = (int)par[0];
	double Pr[2];
  Pr[0] = (double)par[1];
  Pr[1] = (double)par[2];
  double Pr_exit[2];
	Pr_exit[0] = (double)par[3];
	Pr_exit[1] = (double)par[4];
  double Er[2];
	Er[0] = (double)par[5];
	Er[1] = (double)par[6];
  double G0[2], G1[2], G2[2];
	G0[0] = (double)par[7];
	G0[1] = (double)par[8];
	G1[0] = (double)par[9];
	G1[1] = (double)par[10];
	G2[0] = (double)par[11];
	G2[1] = (double)par[12];
	double sign = (double)par[13];
	double Temp = (double)par[14];

  double Scale[3];
	double RatePart[3];    // Stores the rate part for each resonance and the interference term
	double SFactorPart[3]; // Stores the S-factor for each term
	double delta[3];
	
  // double mue = M0*M1/(M0+M1);
  // double R = R0*(pow(M0,1./3.) + pow(M1,1./3.));
  double PEK = 6.56618216E-1 / mue; // a correction factor

	double singular_point;

	// Do the individual resonances
	for(int iRes=0; iRes<2; iRes++){
		double P = PenFactor(x, L[iRes][0], M0, M1, Z0, Z1, R);
		double P_exit, E_exit=0.0;
		double omega = (2. * Jr[iRes] + 1.) / ((2. * J0 + 1.) * (2. * J1 + 1.));

	// TODO Got here!
	
		if (Er[iRes] > 0.0) {
			Scale[0] = P / Pr[iRes];
		} else {
		  Pr[iRes]=0.0;
			Scale[0] = 2.0 * P * 41.80161396 / (mue * gsl_pow_2(R));
		}

		// Exit channel scale
		if(G1[iRes] > 0.0){
			if (Reac.getGamma_index() == 1){
				//std::cout << Exf[iRes] << " " << L[iRes][1] << std::endl;
				Scale[1] = pow((Reac.Q + x - Exf[iRes]) / (Reac.Q + Er[iRes] - Exf[iRes]), (2. * L[iRes][1] + 1.0));
			} else {
				E_exit = Reac.Q + x - Reac.Qexit - Exf[iRes];
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[iRes][1], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[1] = P_exit / Pr_exit[iRes];
				}
			}
		} else {
			Scale[1] = 1.0;
		}

				// Spectator scale (no excitation energy in final state)
		if(G2[iRes] > 0.0){
			if (Reac.getGamma_index() == 2){
				Scale[2] = pow((Reac.Q + x) / (Reac.Q + Er[iRes]), (2. * L[iRes][2] + 1.0));
			} else {
				E_exit = Reac.Q + x - Reac.Qexit;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[iRes][2], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[2] = P_exit / Pr_exit[iRes];
				}
			}
		} else {
			Scale[2] = 1.0;
		}

		/*		// TODO from here!
		for (int i = 1; i < 3; i++) {
			if (G1[iRes] > 0.0) {
				if (i == Reac.getGamma_index()) {
					if (i == 1)
						Scale[i] = pow((Reac.Q + x - Exf[iRes]) /
													 (Reac.Q + Er[iRes] - Exf[iRes]), (2. * L[iRes][i] + 1.0));
					if (i == 2)
						Scale[i] = pow((Reac.Q + x) / (Reac.Q + Er[iRes]), (2. * L[iRes][i] + 1.0));
				} else {
					if (i == 1)
						E_exit = Reac.Q + x - Reac.Qexit - Exf[iRes];
					if (i == 2)
						E_exit = Reac.Q + x - Reac.Qexit;
					if (E_exit > 0.0) {
						P_exit = PenFactor(E_exit, L[iRes][i], M0 + M1 - M2, M2,
															 Z0 + Z1 - Z2, Z2, R);
						Scale[i] = P_exit / Pr_exit[iRes];
					} else {
						Scale[i] = 0.0;
					}
				}
			} else {
				Scale[i] = 1;
			}
		}
		*/

		/*
		std::cout << iRes << "   " << x << " Scale[0]=" << Scale[0] << " Scale[1]=" << Scale[1]
							<< " Scale[2]=" << Scale[2] << std::endl;
		std::cout << "       " << " G0=" << G0[iRes] << " G1=" << G1[iRes]
							<< " G2=" << G2[iRes] << std::endl;
		*/
		
		double S1 = PEK * omega * Scale[0] * G0[iRes] * Scale[1] * G1[iRes];
		double S2 = gsl_pow_2(Er[iRes] - x) +
			0.25 * gsl_pow_2(G0[iRes] * Scale[0] + G1[iRes] * Scale[1] + G2[iRes] * Scale[2]);
		double S3 = exp(-11.605 * x / Temp);

//		std::cout << "            " << "S1=" << S1 << " S2=" << S2 << "S3=" << S3 << std::endl;


		double integrand = S1 * S3 / S2; //*3.7318e10*(pow(mue,-0.5)*pow(Temp,-1.5));
		
		double sfactor = (S1/S2)*exp(0.989534*Z0*Z1*sqrt(mue/x));

		RatePart[iRes] = integrand;
		SFactorPart[iRes] = sfactor;
    delta[iRes] = atan( (G0[iRes]*Scale[0] + G1[iRes]*Scale[1] +
			  G2[iRes]*Scale[2])/(2.0*(x-Er[iRes])));
    if(delta[iRes] < 0.0)delta[iRes]+=M_PI;

		// Check for singular points
		double diff = fabs(x - Er[iRes]);
		singular_point = (diff <= (std::numeric_limits<double>::epsilon() * Er[iRes]));

	}

	// Now the interference term
  // Then, the interference term
  double del = delta[0]-delta[1];
  RatePart[2] = sign*2.0*sqrt(RatePart[0]*RatePart[1])*cos(del);
	SFactorPart[2] = sign*2.0*sqrt(SFactorPart[0]*SFactorPart[1])*cos(del);

  if(isnan(RatePart[2])){
    RatePart[2] = 0.0;
		//    IntfNAN=true;
  }

	double totalIntegrand = RatePart[0] + RatePart[1] - RatePart[2];
	double totalSFactor = SFactorPart[0] + SFactorPart[1] - SFactorPart[2];
	
  if (singular_point)
    totalIntegrand = std::numeric_limits<double>::quiet_NaN();

  // Write the integrand to a file if requested
	//	std::cout << writeIntegrand << " ";
	if (writeIntegrand){
		integrandfile << std::scientific << std::setprecision(9) << x << " " << totalIntegrand
									<< " " << totalSFactor << std::endl;
	}
	

  return totalIntegrand;
}

int Interference::rhs (double x, const double y[], double dydx[], void *params){


  double *par = (double *)params;
  //int writeIntegrand = (int)par[0];
	double Pr[2];
  Pr[0] = (double)par[1];
  Pr[1] = (double)par[2];
  double Pr_exit[2];
	Pr_exit[0] = (double)par[3];
	Pr_exit[1] = (double)par[4];
  double Er[2];
	Er[0] = (double)par[5];
	Er[1] = (double)par[6];
  double G0[2], G1[2], G2[2];
	G0[0] = (double)par[7];
	G0[1] = (double)par[8];
	G1[0] = (double)par[9];
	G1[1] = (double)par[10];
	G2[0] = (double)par[11];
	G2[1] = (double)par[12];
	double sign = (double)par[13];
	double Temp = (double)par[14];


  double Scale[3];
	double RatePart[3];    // Stores the rate part for each resonance and the interference term
	double SFactorPart[3]; // Stores the S-factor for each term
	double delta[3];

  // double mue = M0*M1/(M0+M1);
  // double R = R0*(pow(M0,1./3.) + pow(M1,1./3.));
  double PEK = 6.56618216E-1 / mue; // a correction factor

	//double singular_point;

	// Do the individual resonances
	for(int iRes=0; iRes<2; iRes++){
		double P = PenFactor(x, L[iRes][0], M0, M1, Z0, Z1, R);
		double P_exit, E_exit=0.0;
		double omega = (2. * Jr[iRes] + 1.) / ((2. * J0 + 1.) * (2. * J1 + 1.));

	// TODO Got here!
	
		if (Er[iRes] > 0.0) {
			Scale[0] = P / Pr[iRes];
		} else {
			Scale[0] = 1.0;
			G0[iRes] = 2.0 * P * G0[iRes] * 41.80161396 / (mue * gsl_pow_2(R));
		}

		
		// Exit channel scale
		if(G1[iRes] > 0.0){
			if (Reac.getGamma_index() == 1){
				//std::cout << Exf << " " << L[iRes][1] << std::endl;
				Scale[1] = pow((Reac.Q + x - Exf[iRes]) / (Reac.Q + Er[iRes] - Exf[iRes]), (2. * L[iRes][1] + 1.0));
			} else {
				E_exit = Reac.Q + x - Reac.Qexit - Exf[iRes];
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[iRes][1], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[1] = P_exit / Pr_exit[iRes];
				}
			}
		} else {
			Scale[1] = 1.0;
		}

				// Spectator scale (no excitation energy in final state)
		if(G2[iRes] > 0.0){
			if (Reac.getGamma_index() == 2){
				Scale[2] = pow((Reac.Q + x) / (Reac.Q + Er[iRes]), (2. * L[iRes][2] + 1.0));
			} else {
				E_exit = Reac.Q + x - Reac.Qexit;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[iRes][2], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[2] = P_exit / Pr_exit[iRes];
				}
			}
		} else {
			Scale[2] = 1.0;
		}

		/*
		// TODO from here!
		for (int i = 1; i < 3; i++) {
			if (G1[iRes] > 0.0) {
				if (i == Reac.getGamma_index()) {
					if (i == 1)
						Scale[i] = pow((Reac.Q + x - Exf[iRes]) /
													 (Reac.Q + Er[iRes] - Exf[iRes]), (2. * L[iRes][i] + 1.0));
					if (i == 2)
						Scale[i] = pow((Reac.Q + x) / (Reac.Q + Er[iRes]), (2. * L[iRes][i] + 1.0));
				} else {
					if (i == 1)
						E_exit = Reac.Q + x - Reac.Qexit - Exf[iRes];
					if (i == 2)
						E_exit = Reac.Q + x - Reac.Qexit;
					if (E_exit > 0.0) {
						P_exit = PenFactor(E_exit, L[iRes][i], M0 + M1 - M2, M2,
															 Z0 + Z1 - Z2, Z2, R);
						Scale[i] = P_exit / Pr_exit[iRes];
					} else {
						Scale[i] = 0.0;
					}
				}
			} else {
				Scale[i] = 1;
			}
		}
		*/
		
		double S1 = PEK * omega * Scale[0] * G0[iRes] * Scale[1] * G1[iRes];
		double S2 = gsl_pow_2(Er[iRes] - x) +
			0.25 * gsl_pow_2(G0[iRes] * Scale[0] + G1[iRes] * Scale[1] + G2[iRes] * Scale[2]);
		double S3 = exp(-11.605 * x / Temp);

		double integrand = S1 * S3 / S2; //*3.7318e10*(pow(mue,-0.5)*pow(Temp,-1.5));
		
		double sfactor = (S1/S2)*exp(0.989534*Z0*Z1*sqrt(mue/x));

		RatePart[iRes] = integrand;
		SFactorPart[iRes] = sfactor;
    delta[iRes] = atan( (G0[iRes]*Scale[0] + G1[iRes]*Scale[1] +
			  G2[iRes]*Scale[2])/(2.0*(x-Er[iRes])));
    if(delta[iRes] < 0.0)delta[iRes]+=M_PI;

		// Check for singular points
		//double diff = fabs(x - Er[iRes]);
		//singular_point = (diff <= (std::numeric_limits<double>::epsilon() * Er[iRes]));

	}

	// Now the interference term
  // Then, the interference term
  double del = delta[0]-delta[1];
  RatePart[2] = sign*2.0*sqrt(RatePart[0]*RatePart[1])*cos(del);
	SFactorPart[2] = sign*2.0*sqrt(SFactorPart[0]*SFactorPart[1])*cos(del);

  if(isnan(RatePart[2])){
    RatePart[2] = 0.0;
		//    IntfNAN=true;
  }

	double totalIntegrand = RatePart[0] + RatePart[1] - RatePart[2];
	//	double totalSFactor = SFactorPart[0] + SFactorPart[1] - SFactorPart[2];

	dydx[0] = totalIntegrand;
	
	return GSL_SUCCESS;
}
//}

//----------------------------------------------------------------------
double Interference::getSFactor(double E, int sign, int samp){

  //double PEK = 6.56618216E-1 / mue; // a correction factor

	
	
	double SFactorPart[3], delta[3];
	
	for(int iRes=0; iRes<2; iRes++){

		double E_cm_samp;
		double G0; 
		double G1; 
		double G2;
		double Exf;
		if(samp == -1){
			E_cm_samp = E_cm[iRes];
			G0 = Res[iRes]->getG(0);
			G1 = Res[iRes]->getG(1);
			G2 = Res[iRes]->getG(2);
			Exf = Res[iRes]->getExf();
		} else {
			E_cm_samp = Res[iRes]->getESample(samp);
			G0 = Res[iRes]->getGSample(0,samp);
			G1 = Res[iRes]->getGSample(1,samp);
			G2 = Res[iRes]->getGSample(2,samp);
			sign = sign_sample[samp];
		}
		double P = PenFactor(E, L[iRes][0], M0, M1, Z0, Z1, R);
		double Pr, Pr_exit;
		double P_exit, E_exit = 0.;
		double omega = (2. * Jr[iRes] + 1.) / ((2. * J0 + 1.) * (2. * J1 + 1.));

		double Scale[3];

		double eta = 0.989510*Z0*Z1*sqrt(mue/E);

		/*
			std::cout << iRes << " " << E << " " << G0 << " "
							<< G1 << " " << G2 << "P=" << P << std::endl;
		std::cout << "E_cm_samp=" << E_cm_samp << " Reac.Q=" << Reac.Q
							<< "Reac.Qexit=" << Reac.Qexit << std::endl;
		*/
		
		//  The penetration factor at the resonance energy (the "true" PF)
		if (E_cm_samp > 0.0) {
			Pr = PenFactor(E_cm_samp, L[iRes][0], M0, M1, Z0, Z1, R);
			Scale[0] = P / Pr;
		} else {
			//Scale[0] = 1.0;
			Pr = 0.0;
			Scale[0] = 2.0 * P * 41.80161396 / (mue * gsl_pow_2(R));
		}
		
		// Calculate the exit particle energy, depends on if it is spectator
		if (Reac.getGamma_index() == 2) {
			// if exit particle is observed decay, take final excitation into account
			Pr_exit = PenFactor(E_cm_samp + Reac.Q - Reac.Qexit - Exf,
													L[iRes][1], M0 + M1 - M2, M2,
													Z0 + Z1 - Z2, Z2, R);
		} else if (Reac.getGamma_index() == 1 && NChannels[iRes] == 3) {
			// ignore spectator final excitation if it is spectator
			//			std::cout << E << " E_cm_samp + Reac.Q - Reac.Qexit = " << E_cm_samp + Reac.Q - Reac.Qexit
			//								<< std::endl;
			Pr_exit = PenFactor(E_cm_samp + Reac.Q - Reac.Qexit, L[iRes][2], M0 + M1 - M2, M2,
													Z0 + Z1 - Z2, Z2, R);
		}
		//std::cout << iRes << "   " << E << " Pr=" << Pr << " Pr_exit=" << Pr_exit << std::endl;
		//std::cout << iRes << "   " << E << " Scale[0]=" << Scale[0] << std::endl;
		//		std::cout << Reac.getGamma_index() << std::endl;
		// Exit scale
		if(G1 > 0.0){
			if (Reac.getGamma_index() == 1){
				//std::cout << Exf << " " << L[iRes][1] << std::endl;
				Scale[1] = pow((Reac.Q + E - Exf) / (Reac.Q + E_cm_samp - Exf), (2. * L[iRes][1] + 1.0));
			} else {
				E_exit = Reac.Q + E - Reac.Qexit - Exf;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[iRes][1], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[1] = P_exit / Pr_exit;
				}
			}
		} else {
			Scale[1] = 1.0;
		}
		
		// Spectator scale (no excitation energy in final state)
		if(G2 > 0.0){
			if (Reac.getGamma_index() == 2){
				Scale[2] = pow((Reac.Q + E) / (Reac.Q + E_cm_samp), (2. * L[iRes][2] + 1.0));
			} else {
				E_exit = Reac.Q + E - Reac.Qexit;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[iRes][2], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[2] = P_exit / Pr_exit;
				}
			}
		} else {
			Scale[2] = 1.0;
		}

		/*
		// Exit and spectator scales
		for (int i = 1; i < 3; i++) {
			if (G[iRes][i] > 0.0) {
				if (i == Reac.getGamma_index()) {
					if (i == 1)
						Scale[i] =
              pow((Reac.Q + E - Exf[iRes]) / (Reac.Q + E_cm_samp - Exf[iRes]), (2. * L[iRes][i] + 1.0));
					if (i == 2)
						Scale[i] = pow((Reac.Q + E) / (Reac.Q + E_cm_samp), (2. * L[iRes][i] + 1.0));
				} else {
					if (i == 1)
						E_exit = Reac.Q + E - Reac.Qexit - Exf[iRes];
					if (i == 2)
						E_exit = Reac.Q + E - Reac.Qexit;
					if (E_exit > 0.0) {
						P_exit =
              PenFactor(E_exit, L[iRes][i], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
						Scale[i] = P_exit / Pr_exit;
					} else {
						Scale[i] = 0.0;
					}
				}
			} else {
				Scale[i] = 1;
			}
		}
		*/
		double S1 = exp(eta);
		double S2 = omega * Scale[0] * G0 * Scale[1] * G1;
		double S3 = gsl_pow_2(E_cm_samp - E) +
			0.25 * gsl_pow_2(G0 * Scale[0] + G1 * Scale[1] + G2 * Scale[2]);

		/*
		std::cout << iRes << "   " << E << " Scale[0]=" << Scale[0] << " Scale[1]=" << Scale[1]
							<< " Scale[2]=" << Scale[2] << std::endl;
		std::cout << "       " << " G0=" << G0 << " G1=" << G1
							<< " G3=" << G2 << std::endl;
		*/
		//		std::cout << Scale[0] << " " << Scale[1] << " " << Scale[2] << "\n";
		//std::cout << S1 << " " << S2 << " " << S3 << "\n";
	
		SFactorPart[iRes] = S1 * S2 / S3; //*3.7318e10*(pow(mue,-0.5)*pow(Temp,-1.5));
		delta[iRes] = atan( (G0*Scale[0] + G1*Scale[1] +
												 G2*Scale[2])/(2.0*(E - E_cm_samp)));
		if(delta[iRes] < 0.0)delta[iRes]+=M_PI;
	}

	// Now the interference term
  // Then, the interference term
  double del = delta[0]-delta[1];
	SFactorPart[2] = sign*2.0*sqrt(SFactorPart[0]*SFactorPart[1])*cos(del);
	double Consts = 0.6566/mue;   // pi*hbar^2/(2*mu) in MeV.b
	
	return Consts*(SFactorPart[0] + SFactorPart[1] - SFactorPart[2]);
	
}



void Interference::print() {

  
	std::cout << " Interfering pair " << std::setw(3) << index << std::endl;
	switch(IntfSign){
	case(0):
		std::cout << "  with unknown sign\n";
		break;
	case(1):
		std::cout << "  with positive interference\n";
		break;
	case(-1):
		std::cout << "  with negative interference\n";
		break;
	default:
		std::cout << "  SOME ERROR IN INTERFERENCE!\n";
	}
	std::cout << "  First Resonance:\n";
	Res[0]->print();
	std::cout << "  Second Resonance:\n";
	Res[1]->print();

}

void Interference::write() {

  //  logfile << "--------------------------------------------------" << "\n";
  // logfile << "     This is resonance: " << index << "\n";
  logfile << " Interfering Pair " << std::setw(3) << index << "\n";
	switch(IntfSign){
	case(0):
		logfile << "  with unknown sign\n";
		break;
	case(1):
		logfile << "  with positive interference\n";
		break;
	case(-1):
		logfile << "  with negative interference\n";
		break;
	default:
		logfile << "  SOME ERROR IN INTERFERENCE!\n";
	}
	logfile << "  First Resonance:\n";
	Res[0]->write();
	logfile << "  Second Resonance:\n";
	Res[1]->write();
        /*
	"    E_cm = " << E_cm
    << " +/- " << dE_cm << " MeV\n";
logfile << "                 wg   = " << wg << " +/- " << dwg << " MeV\n";
logfile << "                 Jr   = " << Jr << "\n";
logfile << "                 G1   = " << G[0] << " +/- " << dG[0]
    << " MeV (L = " << L[0] << ")\n";
if (isUpperLimit)
logfile << "                 PT   = " << PT[0] << " +/- " << dPT[0] << "\n";
logfile << "                 G2   = " << G[1] << " +/- " << dG[1]
    << " MeV (L = " << L[1] << ")\n";
if (isUpperLimit)
logfile << "                 PT   = " << PT[1] << " +/- " << dPT[1] << "\n";
logfile << "                 G3   = " << G[2] << " +/- " << dG[2]
    << " MeV (L = " << L[2] << ")\n";
if (isUpperLimit)
logfile << "                 PT   = " << PT[2] << " +/- " << dPT[2] << "\n";
  logfile << "            NChannels = " << NChannels << "\n";
logfile << "                  Exf = " << Exf << "\n";
logfile << "           Integrated = " << isBroad << "\n";
logfile << "          Upper Limit = " << isUpperLimit << "\n";
  logfile << "    Energy Correlated = " << isECorrelated << "\n";
  logfile << "     Width Correlated = " << isWidthCorrelated << "\n";
  logfile << "   Corresponding res. = " << CorresRes << "\n";
  logfile << "                 Frac = " << Frac << "\n";

int NPrintSamples = 5;
logfile << "First " << NPrintSamples << " samples    -------\n";
logfile << "E_cm: ";
//  logfile << E_sample.size() << "\n";
// Print energy samples
for (int s = 0; s < NPrintSamples; s++) {
logfile << E_sample[s] << " ";
}
logfile << "\n";

// Print wg samples
if (wg_sample.size() > 0) {
logfile << "wg: ";
for (int s = 0; s < NPrintSamples; s++) {
logfile << wg_sample[s] << " ";
}
logfile << "\n";
}

// Print Gamma samples
for (int i = 0; i < NChannels; i++) {
if (G_sample[i].size() > 0) {
logfile << "G" << i << ": ";
for (int s = 0; s < NPrintSamples; s++) {
  logfile << G_sample[i][s] << " ";
}
logfile << "\n";
}
}
	*/
  logfile << "\n";
}
