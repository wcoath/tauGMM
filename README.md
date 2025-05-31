# tauGMM
R scripts for running Gaussian mixture models and calculating cutpoints/positivity rates with 95% confidence intervals.

Currently set up for ADNI dataset (UCBERKELEY_TAU_6MM_31May2025.csv) but can be adapted to different datasets easily.
Cutpoints are set to the 99th percentile of the lower (negative) Gaussian. 
These functions can also be applied to similar biomarker datasets that may contain a mix of individuals with and without pathology e.g. amyloid PET.

## installation
Tested on R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.5

R packages used:
tidyverse
mclust
ggpubr
broom

Tested on versions:
ggpubr_0.6.0
mclust_6.1
tidyverse_2.0.0
broom_1.0.6

## set up
Create a sub directory inside tauGMM called raw_data containing an input csv file.

## Output
An 'output' sub-directory should be created by the script with the following:
csv files for cutpoints, positivity rates and individual tau status
pdf plots of bootstrapped cutpoint distributions and Gaussian mixture model fits
