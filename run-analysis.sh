#!/bin/bash
# Run script to run full bacterial CUE analysis
# and produce results figures
# Tom Smith 2021


echo "Removing low quaility data"

Rscript Code/Filtering.R

echo "Fitting Sharpe-Schoolfield model to filtered data"

Rscript Code/SchoolfieldModelling.R

echo "Fitting and bootstrapping BA model to OTR based on Schoolfield fits (bootstrapping will be slow)"

Rscript Code/BAModelling.R

echo "Analysing bootstrap results to investiage CUE TPC directionality"

Rscript Code/CUE-directionality.R

echo "Create the individual plots for figure 1"

Rscript Code/fig1-plots.R

echo "Do the final anyalysis and create the remaining plots from the manuscript"

Rscript Code/make-figures.R

echo "... Finished performing the analyses, output written to /Results/"
