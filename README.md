# RatesMC

This code is the new, improved version of the Monte Carlo reaction
rate code first presented in 2010:

  * Charged-particle thermonuclear reaction rates: I. Monte Carlo
    method and statistical distributions, Longland, R.; Iliadis,
    C.; Champagne, A. E.; Newton, J. R.; Ugalde, C.; Coc, A.;
    Fitzgerald, R.,
    [link](https://ui.adsabs.harvard.edu/abs/2010NuPhA.841....1L/abstract)
  * Charged-particle thermonuclear reaction rates: II. Tables and
    graphs of reaction rates and probability density functions,
    Iliadis, C.; Longland, R.; Champagne, A. E.; Coc, A.;
    Fitzgerald, R.,
    [link](https://ui.adsabs.harvard.edu/abs/2010NuPhA.841...31I/abstract)
  * Charged-particle thermonuclear reaction rates: III. Nuclear
    physics input, Iliadis, C.; Longland, R.; Champagne, A. E.;
    Coc, A.,
    [link](https://ui.adsabs.harvard.edu/abs/2010NuPhA.841..251I/abstract)
  * Charged-particle thermonuclear reaction rates: IV. Comparison to
    previous work, Iliadis, C.; Longland, R.; Champagne, A. E.;
    Coc, A.,
    [link](https://ui.adsabs.harvard.edu/abs/2010NuPhA.841..323I/abstract)
    
Subsequent updates to the method are documented in:

  * Thermonuclear reaction rate of 18Ne(α ,p ) 21Na from Monte Carlo
    calculations, Mohr, P.; Longland, R.; Iliadis, C.,
    [link](https://ui.adsabs.harvard.edu/abs/2014PhRvC..90f5806M/abstract)
  * Correlated uncertainties in Monte Carlo reaction rate
    calculations, Longland, Richard,
    [link](https://ui.adsabs.harvard.edu/abs/2017A%26A...604A..34L/abstract)
  * Correlated energy uncertainties in reaction rate calculations,
    Longland, Richard; de Séréville, Nicolas, 
    [link](https://ui.adsabs.harvard.edu/abs/2020A%26A...642A..41L/abstract)

# Prerequisites

  * GSL - Gnu Scientific Library
  * c++ compiler
  * cmake
  * (openmp - not used at present)

# Installation Instructions

1. Download the code:  
    `git clone https://github.com/rlongland/RatesMC.git`
2. Create a build directory:  
   `cd RatesMC`  
   `mkdir build`  
   `cd build`
3. Compile the code  
   `cmake ..`  
   `make`
   
The code should now be ready to run with  
`./RatesMC`

Update at any time with  
    `git pull`  
    `cd build`  
    `cmake ..`  
    `make`

# Running the code

When you run `./RatesMC`, it will read the input file: `RatesMC.in`,
perform a Monte Carlo reaction rate calculation, and will write
several output files. Some of these are then used by the analysis
scripts below:

* **RatesMC.log**  
  The log file. *Check this for errors*
* **RatesMC.out**  
  The simplified reaction rate output file containing all essential
  rate information
* **RatesMC.full**  
  The Full RatesMC output table. Includes 2-sigma rate uncertainties,
  1-sigma rate uncertainties, Classical Rate,  Median Rate, Mean Rate,
  Log-Normal mu and sigma parameters, and the Anderson-Darling Statistic
* **RatesMC.latex**  
  Latex file for pasting rate into publications
* **RatesMC.cont**  
  Contributions of each resonance to the reaction rate at each temperature
* **RatesMC.integ**  
  Integrals for broad resonances
* **RatesMC.sfact**  
  Astrophysical S-factor for broad resonances and analytical rate
* **RatesMC.samp**  
  Reaction rate samples
* **ParameterSamples.dat**  
  All samples of the input parameters
* **RatesMC.PT**  
  Porter-Thomas samples. Use for diagnostics only

# New features
* **Non-resonant Contributions**  
  Previously non-resonance contributions could be added to the
  reaction rate calculation by including them as a parameterized list
  of `S(0)`, `S'(0)`, and `S''(0)` as explained in Eqn. 3.102 of
  Iliadis' Nuclear Physics of Stars (2nd ed.). However, that requires
  a fit to the astrophysical S-factor. Now, I've included the option
  to include the S-factor directly as a table.   
  There is space for two parts. Each part contains 3 lines:
  1. Energy (in keV)
  2. S (keV.b)
  3. Fractional uncertainty of S-factor (enter a negative value for
  factor uncertainty)  
  For example:  
  ```
  Non-resonant table  
  Enter comma-separated energies (keV) on 1st line, S-factors (keV.b) on 2nd line, and fractional uncertainties on 3rd line
  0      2      4      6      8      10     20     40     60     80     100    200    400    600    800    1000   1500   2000   2500
  10.92  10.81  10.69  10.59  10.48  10.38  9.92   9.15   8.56   8.07   7.67   6.39   5.21   4.65   4.35   4.13   3.76   3.56   3.42
  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04  -1.04
  0  0
  0  0
  0  0
  ```
* **Interfering Resonances**  
  Interferences now work like in the old version of RatesMC. Be warned
  that this feature is still under testing and there is not a lot of
  error checking at present.
* **Correlations**  
  Turning correlations between resonances was a bit of a hidden
  feature before, and was a little finicky with little 'c's to flag
  which resonances should be correlated. This how now been improved:
  at the end of a resonance input, enter 'e' to indicate that it's
  energy-correlated, and 'w' if it's width correlated (including
  resonance strength).   
  Example 1 (resonances at 560 and 877 keV have correlated widths):
  ```
  Ecm    DEcm     wg       Dwg    Jr      G1      DG1     L1      G2      DG2      L2   G3  DG3  L3 Exf   Int
  -43.9  1.5      0        0      1.5     2.2     1.6     1       1.14    0.04      1    0     0    0   0.0   1
  560.0  1.7      0.331    0.041  1.5     6       3       2       0.350   0.056     1    0     0    0   0.0   0  w
  877.0  2.0      3.65     0.51   2.5     11.9    3.6     2       3.1     1.3       1    0     0    0   0.0   0  w

  ```
  Example 2 (resonances at 2101 and 2236 keV have correlated energies):
  ```
  2101.0 50.0     0        0      1.5     210e3  40e3     2       17      8.5       1    48.0  9.2  0   0.0   1  e
  2236.0 50.0     0        0      2.5     240e3  50e3     2       1       0.5       1    520.0 100  0   0.0   1  e
  2426.0 50.0     0        0      3.5     120e3  20e3     2       1       0.5       1    760   140  0   0.0   1
  ```
* **AME masses**  
  Do you want to just read the AME 2020 mass table? Me too! Now,
  insead of entering the particle mass in lines 6-8, you can just type
  the nuclide's name. For example:

  ```
  35Ar(p,g)36K
  ****************************************************************************************************************
  1               ! Zproj
  18              ! Ztarget
  0               ! Zexitparticle (=0 when only 2 channels open)
  1H              ! Aproj   (input either mass in amu, or the nuclide's name, e.g. 23Na)
  35Ar            ! Atarget (input either mass in amu, or the nuclide's name, e.g. 23Na)
  0               ! Aexitparticle (=0 when only 2 channels open)
  ```
  
* **AME atomic number**  
  Keep forgetting the atomic number of your nuclei? Enter the
  nuclide's name in the charge section. Yes I know the isotope is
  irrelevent.

  ```
  35Ar(p,g)36K
  ****************************************************************************************************************
  1H              ! Zproj
  35Ar            ! Ztarget
  0               ! Zexitparticle (=0 when only 2 channels open)
  1H              ! Aproj   (input either mass in amu, or the nuclide's name, e.g. 23Na)
  35Ar            ! Atarget (input either mass in amu, or the nuclide's name, e.g. 23Na)
  0               ! Aexitparticle (=0 when only 2 channels open)
  ```

* **Nubase ground-state spin**
  Tired of looking up the ground state spin of your reacting nuclides?
  Guess what, you can enter their names to use the NuBase 2020
  evaluation. This *ignores parentheses* so use with caution!
  
  ```
  35Ar(p,g)36K
  ****************************************************************************************************************
  1               ! Zproj     (or nuclide's name to use AME 2020)
  18              ! Ztarget
  0               ! Zexitparticle (=0 when only 2 channels open)
  1.00728         ! Aproj     (input either nuclear mass in amu, or the nuclide's name, e.g. 23Na)
  34.9654         ! Atarget
  0               ! Aexitparticle (=0 when only 2 channels open)
  1H              ! Jproj     (or nuclide's name to use NuBase 2020)
  35Ar            ! Jtarget
  ```
  
* **Separation energies** 
  Sick of calculating separation energies from
  nuclide masses? Yep, you guessed it, RatesMC can finally do that for
  you! 
  ```
  35Ar(p,g)36K
  ****************************************************************************************************************
  1H              ! Zproj     (or nuclide's name to use AME 2020)
  35Ar            ! Ztarget
  0               ! Zexitparticle (=0 when only 2 channels open)
  1H              ! Aproj     (input either nuclear mass in amu, or the nuclide's name, e.g. 23Na)
  35Ar            ! Atarget
  0               ! Aexitparticle (=0 when only 2 channels open)
  1H              ! Jproj     (or nuclide's name to use NuBase 2020)
  35Ar            ! Jtarget
  0.0             ! Jexitparticle (=0 when only 2 channels open)
  AME             ! projectile separation energy (keV) - Enter 'AME' for auto calculation
  AME             ! exit particle separation energy (=0 when only 2 channels open)

  ```
* **Factor uncertainties**  
  Factor uncertainties on parameters can finally be input directly
  into the code! **Beware: this only works for partial widths and
  Porter Thomas mean values, not energies**.  
  For example, if you want to enter a factor of two uncertainty for
  the first partial width of a resonance, you would enter the
  following:
  ```
  Ecm     DEcm    wg      Dwg     J     G1        DG1        L1    G2      DG2      L2   G3  DG3  L3  Exf   Int
  200     2       0       0       2     0.0372    -2.0       1     2.7e-4  1.35e-4  1    0   0    0   0.0   1
  ```
# Analysis scripts
Most of these will need to be edited to make sure the correct files
are being analysed.

* **PlotUncertainties.R**  
  New code to plot the 1, 2, and 3-sigma uncertainties of the rate, as
  well as a comparison to literature if desired
* **PlotCompare.R**  
  Plots the 1-sigma uncertainties of the present rate and a literature
  rate
* **PlotContribution.R**  
  Plots the contribution bands of each resonance as a function of
  temperature
* **PlotIntegrand.R**  
  Plots the integrands for broad resonances. Be sure to adjust the
  temperature you're looking for at the top of the code
* **PlotSFactor.R**  
  Plots the astrophysical s-factor for all resonances. This has now
  been updated from the previous version to work well! 
* **PlotPanel6.R**  
  Plots the reaction rate probability density distribution for 6
  useful temperatures
* **PlotPanelall.R**  
  Plots the reaction rate probability density distribution for all
  temperatures
* **PlotParameter.R**  
  Diagnostics. Check the probability distribution of input parameters
* **PlotCorrelations.R**  
  Diagnostics. Show correlations of rate (at a given temperature) with
  all input parameters. Also show corner plot correlations between
  everything (should be mostly uncorrelated unless energy or width
  correlations are turned on)
* **PlotPT.R**  
  Diagnostics. Plot the Porter-Thomas samples
  
# Licence

Copyright (C) 2022  Richard Longland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Author: R. Longland  
Email: rllongla@ncsu.edu
