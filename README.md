# TempDepBactCUE

## Table of Contents
1. Overview
2. Instructions for Use
    1. Software requirements
    2. Full pipeline
3. Repository Contents
    1. Data
    2. Code
    3. Results
4. Citation and license

## Overview

Temperature dependence of bacterial carbon use efficiency

Replicate the analyses for our paper: "Systematic variation in the temperature dependence of bacterial carbon use efficiency", doi: [10.1101/2020.09.14.296095](https://doi.org/10.1101/2020.09.14.296095)

## Instructions for use

### Software requirements

Analysis pipeline was built and tested in R version 3.6.3

### Full pipeline

The whole analysis pipeline to automatically install required R packages, clean and fit models to raw data then analyse the results can be run using the shell script from the main directory `bash run-analysis.sh`

Note that because this code pipeline uses a bootstrapping method, there are likely to be small changes in the numbers reported when you run this through compared to those reported in our manuscript, 
however our results are robust to this and conclusions are expected to be qualitatively unchanged with each run.


## Repository contents

### Data

* dataset_rates_taxonomy.csv - these are the raw rates data
* isolates_growth_temperatures.csv - these are details of isolates used and some assosciated functional traits
* filtered_data.csv - output of Filtering.R; this is the raw data used in TPC fitting
* group_means_non_aggregated_fluxes and all_means_non_aggregated - these are the data synthesis TPCs from Smith et al. 2019 (doi: 10.1038/s41467-019-13109-1), pulled from that project's github repository: https://github.com/smithtp/hotterbetterprokaryotes

### Code

* Packages.R - automatically install and load required R packages
* Filtering.R - take raw data, subset to aerobes, and filter out unreliable rate estimates
* SchoolfieldModelling.R - fit Schoolfield model to filtered data
* BAModelling.R - fit and bootstrap Boltzmann-Arrhenius model to OTR based on Schoolfield fits
* CUE-directionality.R - check bootstrapped CI overlaps to test CUE directionality; write results
* fig1-plots.R - create all of the individual CUE TPC plots for figure 1
* make-figures.R - create Figs 2 and 3 of the manuscript.


### Results (These will be built by the code pipeline)

* Schoolfield_growth.csv - results of Schoolfield fits to growth rate data
* Schoolfield_resp.csv - results of Schoolfield fits to respiration rate data
* Schoolfield_results.csv - combined growth and respiration rate fits, filtered to good estimates for further analysis
* BA_fits_phase1.csv - fits of Boltzmann-Arrhenius to the OTR (determined by the Schoolfield fits)
* bootstrapped_BA_fits.csv - full output of bootstrapping the BA fits
* bootstrap_results_BA.csv - summary of bootstrapped fits
* CUE_directionality.csv - output of checking CI overlaps in CUE-directionality.R, for supp table 2


## Citation and license

This software is licensed under a GNU General Public License v3.0. If using this data or analysis code, please cite our bioRXiv preprint (https://doi.org/10.1101/2020.09.14.296095)
