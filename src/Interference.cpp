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

Interference::Interference(Reaction &R)
	: Reac(R) {}

Interference::~Interference() {}


//----------------------------------------------------------------------
// Function to numerically integrate broad resonances
double Interference::calcBroad(double T) {

  double classicalRate = 0.0; //, ARate;
  /*
  // Calculate the rate samples
  for (int s = 0; s < NSamples; s++) {

  if (s % 10 == 0) {
  // \r goes back to the beginning of the line.
  std::cout << "\r" << 100 * s / NSamples << "% Complete for Resonance "
            << index + 1 << std::flush;
  }

            // Calculate the single integrated rate sample
  Rate_sample[s] = NumericalRate(T, E_sample[s], G_sample[0][s],
                               G_sample[1][s], G_sample[2][s], erFrac[0][s],
                               erFrac[1][s], erFrac[2][s], false);

  }

  // And the central value, which is the classical rate
    // Write the integrand to a file
  classicalRate = NumericalRate(T, E_cm, G[0], G[1], G[2], 1.0, 1.0, 1.0, true);
    integrandfile << std::endl;
          */
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
double Interference::NumericalRate(double T, double E, double G0, double G1,
                                double G2, double erFrac0, double erFrac1,
                                double erFrac2, bool writeIntegrand) {

  double ARate = 0.0;

  double Pr(0.0), Pr_exit(0.0);

  InterferencePtr = this;

  // if reaction is endothermic, need to make minimum energy
  // enough to ensure no integration over negative energies
  double E_min = EMin;

  // if particle is in spectator channel, integration should
  //  not be truncated
  if (Reac.Qexit > Reac.Q && Reac.getGamma_index() == 2)
    E_min += Reac.Qexit + Exf - Reac.Q;

  double E_max = 10.0;

  //  ofstream evsr;
  //  ofstream testhist;
  //  testhist.open("integrands.dat");

	// Scale the partial widths by the energy effect
	G0 *= erFrac0;
	G1 *= erFrac1;
	G2 *= erFrac2;
	
  // IF we input E>0, but the sample is <0, we need to treat it as a
  // subthreshold resonance. To do this, we need to convert
  // G into C2S*Theta_sp
  if (E_cm > 0.0 && E < 0.0) {
    ErrorFlag = true;
    SampledNegCount++;
    G0 = mue * gsl_pow_2(R) * G0 /
         (2.0 * 41.80161396 * PenFactor(E_cm, L[0], M0, M1, Z0, Z1, R));
    //    std::cout << "Positive resonance went negative!\n";
    //    std::cout << "E_cm = " << E_cm << "E_sample = " << E << "G[0] = " <<
    //    G[0]
    //              << " G_sample = " << G0 << std::endl;
  } else if (E_cm < 0.0 && E > 0.0) {
    // or convert to real resonance if E>0.0
    ErrorFlag = true;
    SubSampledPosCount++;
    G0 = G0 * 2.0 * 41.80161396 * PenFactor(E, L[0], M0, M1, Z0, Z1, R) /
         (mue * gsl_pow_2(R));
		//		std::cout << "Negative resonance went positive!\n";
		//		std::cout << "E_cm = " << E_cm << "E_sample = " << E << "G[0] = " << G[0]
		//							<< " G_sample = " << G0 << std::endl;
  }

  //  The penetration factor at the resonance energy (the "true" PF)
  if (E > 0.0) {
    Pr = PenFactor(E, L[0], M0, M1, Z0, Z1, R);
    //    std::cout << "Pr = " << Pr << "\n";
		if(isZero(Pr))return 0.0;
	} else {
    Pr = 0.0;
  }
	//std::cout << "Pr = " << Pr << "\n";

	
  // Calculate the exit particle energy, depends on if it is spectator
  // if(NChannels[j]==3){
  if (Reac.getGamma_index() == 2) {
    // if exit particle is observed decay, take final excitation into account
    Pr_exit = PenFactor(E + Reac.Q - Reac.Qexit - Exf, L[1], M0 + M1 - M2, M2,
                        Z0 + Z1 - Z2, Z2, R);
    // cout << "Exit energy = " << E+Reac.Q-Reac.Qexit-Exf << endl;
  } else if (Reac.getGamma_index() == 1 && NChannels == 3) {
    // ignore spectator final excitation if it is spectator
    Pr_exit = PenFactor(E + Reac.Q - Reac.Qexit, L[2], M0 + M1 - M2, M2,
                        Z0 + Z1 - Z2, Z2, R);
    // cout << "Exit energy = " << E+Reac.Q-Reac.Qexit << endl;
  }

  //--------------------------------------------------
  // GSL Integration functions
  double result, error;
  size_t nevals;

  // Define integration limits
  //double x = E_min, x1 = E_max;
  //double pole = E;

  // The array that needs to be passed to the integration function
  double alpha[8];
  alpha[0] = (int)writeIntegrand;
  alpha[1] = Pr;
  alpha[2] = Pr_exit;
  alpha[3] = E;
  alpha[4] = T;
  alpha[5] = G0;
  alpha[6] = G1;
  alpha[7] = G2;

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

  gsl_integration_cquad_workspace *w =
      gsl_integration_cquad_workspace_alloc(10000);
  int status = gsl_integration_cquad(&F,      // Function to be integrated
                                     E_min,   // Where known singularity is
                                     E_max,   // number of singularities
                                     1e-100,   // absolute error
                                     1e-6,    // relative error
                                     w,       // workspace
                                     &result, // The result
                                     &error, &nevals);
  gsl_integration_cquad_workspace_free(w);

	//status = -1;
  // If the integration errored, use the slower ODE method
  if (status != 0) {
	//	std::cout << "Integration error = " << gsl_strerror(status) << "\n";
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
      if(x < E && (x + h) > (E-(3.0*(G0+G1+G2))))smallstep=true;
			
      // Now make step small if it isn't already, can be fairly large because
      // it starts on the resonance wing
      if(smallstep && h>(G0+G1+G2)){
				h = (G0+G1+G2);
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
  double Pr = (double)par[1];
  double Pr_exit = (double)par[2];
  double Er = (double)par[3];
  double Temp = (double)par[4];
  double G0 = (double)par[5];
  double G1 = (double)par[6];
  double G2 = (double)par[7];

  //  this->print();

  // std::cout << Pr << " " << Pr_exit << " " << Er << " "
  //	    << Temp << " " << G0 << " " << G1 << " " << G2 << " " << "\n";

  // std::cout << "R = " << R << "\n";

  double Scale[3];

  // double mue = M0*M1/(M0+M1);
  // double R = R0*(pow(M0,1./3.) + pow(M1,1./3.));
  double PEK = 6.56618216E-1 / mue; // a correction factor
  double P = PenFactor(x, L[0], M0, M1, Z0, Z1, R);
  double P_exit, E_exit = 0.;
  double omega = (2. * Jr + 1.) / ((2. * J0 + 1.) * (2. * J1 + 1.));

  // cout << Jr << " " << J0 << " " << J1 << " " << mue << " " << R << " " <<
  // PEK << " " << omega << "\n";

  // std::cout << P << " " << omega << "\n";

  if (Er > 0.0) {
    Scale[0] = P / Pr;
  } else {
    Scale[0] = 1.0;
    G0 = 2.0 * P * G0 * 41.80161396 / (mue * gsl_pow_2(R));
  }
  // std::cout << "Scale[0] = " << Scale[0] << " E = " << x << " G0(E) = " << G0
  // << "\n";

  for (int i = 1; i < 3; i++) {
    if (G[i] > 0.0) {
      if (i == Reac.getGamma_index()) {
        if (i == 1)
          Scale[i] =
              pow((Reac.Q + x - Exf) / (Reac.Q + Er - Exf), (2. * L[i] + 1.0));
        if (i == 2)
          Scale[i] = pow((Reac.Q + x) / (Reac.Q + Er), (2. * L[i] + 1.0));
      } else {
        if (i == 1)
          E_exit = Reac.Q + x - Reac.Qexit - Exf;
        if (i == 2)
          E_exit = Reac.Q + x - Reac.Qexit;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[i], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[i] = P_exit / Pr_exit;
        } else {
          Scale[i] = 0.0;
        }
      }
    } else {
      Scale[i] = 1;
    }
  }

  double S1 = PEK * omega * Scale[0] * G0 * Scale[1] * G1;
  double S2 = gsl_pow_2(Er - x) +
              0.25 * gsl_pow_2(G0 * Scale[0] + G1 * Scale[1] + G2 * Scale[2]);
  double S3 = exp(-11.605 * x / Temp);

  double integrand = S1 * S3 / S2; //*3.7318e10*(pow(mue,-0.5)*pow(Temp,-1.5));

	double sfactor = (S1/S2)*exp(0.989534*Z0*Z1*sqrt(mue/x)); 

  //  if(integrand < 1.e-99)integrand=0.0;

  //  std::cout << x << " " << integrand << "\n";

  // cout << x << "\t" << dydx[0] << endl;

  //  integrand = gsl_max(integrand,1e-300);
  double diff = fabs(x - Er);
  bool singular_point = (diff <= (std::numeric_limits<double>::epsilon() * Er));

  if (singular_point)
    integrand = std::numeric_limits<double>::quiet_NaN();

  // Write the integrand to a file if requested
	//	std::cout << writeIntegrand << " ";
	if (writeIntegrand){
		integrandfile << std::scientific << std::setprecision(9) << x << " " << integrand
									<< " " << sfactor << " " << P << std::endl;
	}
	
  // astrohpysical s-factor
  // cout << x << "\t" << x*(S1/S2)/exp(-0.989534*Z0*Z1*sqrt(mue/x)) << endl;
  // Numerator of x section
  // cout << x << "\t" << Scale[0]*G0*Scale[1]*G1 << endl;
  // Denominator of x-section
  // cout << x << "\t" << S2 << "\t" << gsl_pow_2(Er-x) << endl;
  // G0+G1 and G2 and G0+G1+G2
  // cout << x << "\t" << G0*Scale[0]+G1*Scale[1] << "\t" << G2*Scale[2] << "\t"
  // <<
  //  G0*Scale[0]+G1*Scale[1]+G2*Scale[2] << endl;
  // all
  // cout << x << "\t" << x*(S1/S2)/exp(-0.989534*Z0*Z1*sqrt(mue/x)) << "\t" <<
  // S2 <<
  //  "\t" << G0*Scale[0]+G1*Scale[1] << "\t" << G2*Scale[2] << "\t" <<
  //  G0*Scale[0]+G1*Scale[1]+G2*Scale[2] << endl;
  // penetration factor for exit paritlce
  //  cout << x << "\t" << P_exit << endl;
  // std::cout << integrand << "\n";

  return integrand;
}

int Interference::rhs (double x, const double y[], double dydx[], void *params){

	double *par = (double *)params;
  double Pr = (double)par[1];
  double Pr_exit = (double)par[2];
  double Er = (double)par[3];
  double Temp = (double)par[4];
  double G0 = (double)par[5];
  double G1 = (double)par[6];
  double G2 = (double)par[7];

  //  this->print();

  // std::cout << Pr << " " << Pr_exit << " " << Er << " "
  //	    << Temp << " " << G0 << " " << G1 << " " << G2 << " " << "\n";

  // std::cout << "R = " << R << "\n";

  double Scale[3];

  //double mue = M0*M1/(M0+M1);
  // double R = R0*(pow(M0,1./3.) + pow(M1,1./3.));
  double PEK = 6.56618216E-1 / mue; // a correction factor
  double P = PenFactor(x, L[0], M0, M1, Z0, Z1, R);
  double P_exit, E_exit = 0.;
  double omega = (2. * Jr + 1.) / ((2. * J0 + 1.) * (2. * J1 + 1.));

  // cout << Jr << " " << J0 << " " << J1 << " " << mue << " " << R << " " <<
  // PEK << " " << omega << "\n";

  // std::cout << P << " " << omega << "\n";

  if (Er > 0.0) {
    Scale[0] = P / Pr;
  } else {
    Scale[0] = 1.0;
    G0 = 2.0 * P * G0 * 41.80161396 / (mue * gsl_pow_2(R));
  }
  // std::cout << "Scale[0] = " << Scale[0] << " E = " << x << " G0(E) = " << G0
  // << "\n";

  for (int i = 1; i < 3; i++) {
    if (G[i] > 0.0) {
      if (i == Reac.getGamma_index()) {
        if (i == 1)
          Scale[i] =
              pow((Reac.Q + x - Exf) / (Reac.Q + Er - Exf), (2. * L[i] + 1.0));
        if (i == 2)
          Scale[i] = pow((Reac.Q + x) / (Reac.Q + Er), (2. * L[i] + 1.0));
      } else {
        if (i == 1)
          E_exit = Reac.Q + x - Reac.Qexit - Exf;
        if (i == 2)
          E_exit = Reac.Q + x - Reac.Qexit;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[i], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[i] = P_exit / Pr_exit;
        } else {
          Scale[i] = 0.0;
        }
      }
    } else {
      Scale[i] = 1;
    }
  }

  double S1 = PEK * omega * Scale[0] * G0 * Scale[1] * G1;
  double S2 = gsl_pow_2(Er - x) +
              0.25 * gsl_pow_2(G0 * Scale[0] + G1 * Scale[1] + G2 * Scale[2]);
  double S3 = exp(-11.605 * x / Temp);

  double integrand = S1 * S3 / S2; //*3.7318e10*(pow(mue,-0.5)*pow(Temp,-1.5));

	dydx[0] = integrand;
	
	return GSL_SUCCESS;
}
//}

//----------------------------------------------------------------------
double Interference::getSFactor(double E){

  //double PEK = 6.56618216E-1 / mue; // a correction factor
  double P = PenFactor(E, L[0], M0, M1, Z0, Z1, R);
	double Pr, Pr_exit;
  double P_exit, E_exit = 0.;
  double omega = (2. * Jr + 1.) / ((2. * J0 + 1.) * (2. * J1 + 1.));

	double Scale[3];

  double eta = 0.989510*Z0*Z1*sqrt(mue/E);

	
  //  The penetration factor at the resonance energy (the "true" PF)
  if (E_cm > 0.0) {
    Pr = PenFactor(E_cm, L[0], M0, M1, Z0, Z1, R);
  } else {
    Pr = 0.0;
  }

  // Calculate the exit particle energy, depends on if it is spectator
  if (Reac.getGamma_index() == 2) {
    // if exit particle is observed decay, take final excitation into account
    Pr_exit = PenFactor(E_cm + Reac.Q - Reac.Qexit - Exf, L[1], M0 + M1 - M2, M2,
                        Z0 + Z1 - Z2, Z2, R);
  } else if (Reac.getGamma_index() == 1 && NChannels == 3) {
    // ignore spectator final excitation if it is spectator
    Pr_exit = PenFactor(E_cm + Reac.Q - Reac.Qexit, L[2], M0 + M1 - M2, M2,
                        Z0 + Z1 - Z2, Z2, R);
  }

	// Entrance particle scale
  if (E_cm > 0.0) {
    Scale[0] = P / Pr;
  } else {
    Scale[0] = 2.0 * P * 41.80161396 / (mue * gsl_pow_2(R));
		//Scale[0] = 1.0;
    //		G0 = 2.0 * P * G0 * 41.80161396 / (mue * gsl_pow_2(R));
  }

	// Exit and spectator scales
  for (int i = 1; i < 3; i++) {
    if (G[i] > 0.0) {
      if (i == Reac.getGamma_index()) {
        if (i == 1)
          Scale[i] =
              pow((Reac.Q + E - Exf) / (Reac.Q + E_cm - Exf), (2. * L[i] + 1.0));
        if (i == 2)
          Scale[i] = pow((Reac.Q + E) / (Reac.Q + E_cm), (2. * L[i] + 1.0));
      } else {
        if (i == 1)
          E_exit = Reac.Q + E - Reac.Qexit - Exf;
        if (i == 2)
          E_exit = Reac.Q + E - Reac.Qexit;
        if (E_exit > 0.0) {
          P_exit =
              PenFactor(E_exit, L[i], M0 + M1 - M2, M2, Z0 + Z1 - Z2, Z2, R);
          Scale[i] = P_exit / Pr_exit;
        } else {
          Scale[i] = 0.0;
        }
      }
    } else {
      Scale[i] = 1;
    }
  }

	double Consts = 0.6566/mue;   // pi*hbar^2/(2*mu) in MeV.b
	double S1 = exp(eta);
  double S2 = omega * Scale[0] * G[0] * Scale[1] * G[1];
  double S3 = gsl_pow_2(E_cm - E) +
              0.25 * gsl_pow_2(G[0] * Scale[0] + G[1] * Scale[1] + G[2] * Scale[2]);

	//	std::cout << S1 << " " << S2 << " " << S3 << "\n";
	
  double SFactor = Consts * S1 * S2 / S3; //*3.7318e10*(pow(mue,-0.5)*pow(Temp,-1.5));

	return SFactor;
	
}



void Interference::print() {

  
	std::cout << " Interfering pair \n";
	/*
	cout << std::setw(3) << index << "    E_cm = " << E_cm
 << " +/- " << dE_cm << "\n";
cout << "                 wg   = " << wg << " +/- " << dwg << "\n";
cout << "                 Jr   = " << Jr << "\n";
cout << "                 G1   = " << G[0] << " +/- " << dG[0]
 << " (L = " << L[0] << ")\n";
if (isUpperLimit)
cout << "                 PT   = " << PT[0] << " +/- " << dPT[0] << "\n";
cout << "                 G2   = " << G[1] << " +/- " << dG[1]
 << " (L = " << L[1] << ")\n";
if (isUpperLimit)
cout << "                 PT   = " << PT[1] << " +/- " << dPT[1] << "\n";
cout << "                 G3   = " << G[2] << " +/- " << dG[2]
 << " (L = " << L[2] << ")\n";
if (isUpperLimit)
cout << "                 PT   = " << PT[2] << " +/- " << dPT[2] << "\n";
cout << "                 Exf  = " << Exf << "\n";
cout << "           Integrated = " << isBroad << "\n";
cout << "          Upper Limit = " << isUpperLimit << "\n";
  cout << "    Energy Correlated = " << isECorrelated << "\n";
  cout << "     Width Correlated = " << isWidthCorrelated << "\n";
  cout << "   Corresponding res. = " << CorresRes << "\n";
  cout << "                 Frac = " << Frac << "\n";

//  cout << "--------------------------------------------------" << "\n";
int NPrintSamples = 5;
cout << "First " << NPrintSamples << " samples    -------\n";
cout << "E_cm: ";
//  cout << E_sample.size() << "\n";
// Print energy samples
for (int s = 0; s < NPrintSamples; s++) {
cout << E_sample[s] << " ";
}
cout << "\n";

// Print wg samples
if (wg_sample.size() > 0) {
cout << "wg: ";
for (int s = 0; s < NPrintSamples; s++) {
cout << wg_sample[s] << " ";
}
cout << "\n";
}

// Print Gamma samples
for (int i = 0; i < NChannels; i++) {
if (G_sample[i].size() > 0) {
cout << "G" << i << ": ";
for (int s = 0; s < NPrintSamples; s++) {
  cout << G_sample[i][s] << " ";
}
cout << "\n";
}
}
	*/
	std::cout << "\n";
}
void Interference::write() {

  //  logfile << "--------------------------------------------------" << "\n";
  // logfile << "     This is resonance: " << index << "\n";
  logfile << " Interfering Pair \n";
  /*
  logfile << std::setw(3) << index << "    E_cm = " << E_cm
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
