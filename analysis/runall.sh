#!/bin/bash
# A script to batch run RatesMC and main analysis R scripts

## Ask user to make sure TMatch is in place
read -p "Are you using the correct TMatch.HF file? Press [Enter] to continue "

## Run RatesMC
(trap 'kill 0' SIGINT; ./RatesMC & (sleep 5; Rscript PlotSFactor.R) & wait)

## Do the HF matching
Rscript TMatch.R

## Make the LaTeX tables
Rscript makeLaTeX.R

## Contribution plot
Rscript PlotContribution.R

## Uncertainties
Rscript PlotUncertainties.R

## Rate histograms
Rscript PlotPanel6.R
Rscript PlotPanelall.R

## Astrophysical S-Factor
Rscript PlotSFactor.R
