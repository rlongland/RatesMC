#!/bin/bash
# A script to batch run RatesMC and main analysis R scripts

## Run RatesMC
./RatesMC

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
##Rscript PlotSFactor.R

