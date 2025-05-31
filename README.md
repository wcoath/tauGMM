# tauGMM

R scripts for running Gaussian mixture models and calculating cutpoints, positivity status and rates with bootstrapped 95% confidence intervals.

- Currently set up for ADNI dataset (UCBERKELEY_TAU_6MM_31May2025.csv) but can be adapted to different datasets easily.
- Cutpoints are set to the 99th percentile of the lower (negative) Gaussian. This can be altered by changing the zscore in GMM_tau_PET_functions.R 
- These functions can also be applied to similar biomarker datasets that are thought to contain a mix of individuals with and without pathology e.g. amyloid PET.

## Set up/installation

- Clone this repository
- Create a sub directory inside tauGMM called 'raw_data' containing an input csv file
- Update paths and data preparation lines in run_GMM_tau_PET.R

### R set up
- Tested on R version 4.3.2 (2023-10-31).
- Platform: aarch64-apple-darwin20 (64-bit)
- Running under: macOS 15.5

Install required R packages:
- tidyverse (tested with v2.0.0)
- mclust (v6.1)
- ggpubr (v0.6.0)
- broom (v1.0.6)

## Output
- Run the run_GMM_tau_PET.R script
An 'output' sub-directory should be created by the script containing the following:
 - csv files for cutpoints, positivity rates and individual tau status
 - pdf plots of bootstrapped cutpoint distributions and Gaussian mixture model fits
