# Author: Will Coath
# Date: 20250531

# run Gaussian mixture models and calc cutpoints/positivity - set up for ADNI data (UCBERKELEY_TAU_6MM_31May2025.csv)
# cutpoints are set to the 99th percentile of the lower negative Gaussian
# i have a sub-directory inside tauGMM called raw_data containing a tau PET csv file

# this should create an output sub-directory and
# save csv files for cutpoints, positivity rates and individual tau status
# save plots of bootstrapped cutpoints and Gaussian mixture model fits

# libraries ####
library(tidyverse)
library(mclust)
library(ggpubr)
library(broom)

# source functions
source("GMM_tau_PET_functions.R") #update path if not in same dir

# set a seed
s <- 8235235

# define paths ####

# specify root path
root_dir <- file.path("path","to","your","tauGMM")

# create new output sub-directory
out_dir <- file.path(root_dir,"output") 
dir.create(out_dir) #if already created should just print warning

# specify input csv 
infile <- file.path(root_dir,"raw_data","UCBERKELEY_TAU_6MM_31May2025.csv")

# import data and select regions ####
# modify to select whatever regions you would like

df <- read.csv(file = infile) %>%
  filter(TRACER=="FTP",VISCODE=="init" | VISCODE=="bl") %>%
  dplyr::select(LONIUID:TRACER,!contains("_LH_") & !contains("_RH_") & contains("CTX_") & contains("SUVR")) %>%
  drop_na()

# pivot longer by region
mod.df <- df %>%
  pivot_longer(.,cols=CTX_ENTORHINAL_SUVR:CTX_TRANSVERSETEMPORAL_SUVR,names_to = "region",values_to = "suvr") %>%
  mutate(region=gsub("_SUVR","",region))

set.seed(s)

# run with 1:3 Gaussians and see what model is best fit for each region 
# mclust uses BIC to decide on the best fitting region

cutpoints.df <- mod.df %>% 
  dplyr::select(region,suvr) %>% # only want these two columns 
  group_by(region) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ Mclust(.x$suvr,1:3,verbose = FALSE))) %>%
  mutate(s = purrr::map(fit, summary)) %>%
  mutate(params =  purrr::map(s, get_params)) %>%
  unnest_wider(params) %>%
  dplyr::select(-c(data,fit,s)) %>%
  ungroup()

cutpoints.df$n_gaussians <- as.factor(cutpoints.df$n_gaussians)
cutpoints.df$mod_type <- as.factor(cutpoints.df$mod_type)

# now for consistency restrict to 2 Gaussian with unequal variance
# then perform bootstrapping to get cutpoint uncertainty

set.seed(s)
cutpoints2k.df <- mod.df %>% 
  dplyr::select(region,suvr) %>% # only want these two columns 
  group_by(region) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ Mclust(.x$suvr,G = 2, modelNames = "V",verbose = FALSE))) %>%
  mutate(s = purrr::map(fit, summary)) %>%
  mutate(params =  purrr::map(s, get_params)) %>%
  mutate(bootmod = purrr::map(fit, ~ MclustBootstrap(.x, nboot = 5000, verbose = FALSE))) %>%
  mutate(bootparams =  purrr::map(bootmod, get_bootcutpoints)) %>%
  mutate(boot_plot = purrr::map(bootparams, get_bootplot)) %>%
  mutate(bootparams =  purrr::map(bootparams, remove_cutpoints)) %>%
  unnest_wider(params) %>%
  unnest_wider(bootparams) %>%
  dplyr::select(-c(data,fit,s,bootmod)) %>%
  ungroup()

cutpoints2k.df$n_gaussians <- as.factor(cutpoints2k.df$n_gaussians)
cutpoints2k.df$mod_type <- as.factor(cutpoints2k.df$mod_type)
cutpoints2k.df$orig_n_gaussians <- cutpoints.df$n_gaussians # add best fitting model result
cutpoints2k.df$orig_mod_type <- cutpoints.df$mod_type

# calculate positivity status
status2k.df <- left_join(mod.df,cutpoints2k.df) %>%
  dplyr::select(-boot_plot) %>%
  mutate(status = ifelse(suvr >= cutpoint,'Positive','Negative')) %>% # cutpoint estimate
  mutate(status_low = ifelse(suvr >= boot_high,'Positive','Negative')) %>% # lower CI cutpoint, for higher status estimate
  mutate(status_high = ifelse(suvr >= boot_low,'Positive','Negative')) # higher CI cutpoint, for lower status estimate

# calculate positivity rates by region
rates2k.df <- status2k.df %>% 
  dplyr::select(region,status,status_low,status_high) %>%
  group_by(region) %>%
  nest() %>%
  mutate(rates = purrr::map(data,get_rates)) %>%
  unnest_wider(rates) %>%
  dplyr::select(-c(data)) %>%
  ungroup()

# remove GMM parameters from cutpoints and merge with rates
rates2k.df <- cutpoints2k.df %>%
  dplyr::select(-c(mu1:mp3,boot_plot)) %>%
  left_join(.,rates2k.df)

# save bootstrap cutpoint distribution plots ####
for (i in 1:length(cutpoints2k.df$region)) {
  r <- as.character(cutpoints2k.df$region[[i]])
  p <- cutpoints2k.df$boot_plot[[i]]
  print(paste0("saving bootstrapped cutpoint plot: ",r))
  pdf(paste0(out_dir,'/','GMM_bootstrap_cutpoint_plot_',r,'.pdf'), width=6.5, height=5)
  print(p)
  dev.off()
}

# remove bootplots from data frame
cutpoints2k.df <- cutpoints2k.df %>%
  dplyr::select(-boot_plot)

# save mixture model plots ####
for (i in 1:length(cutpoints2k.df$region)) {
  r <- as.character(cutpoints2k.df$region[[i]])
  p <- get_mixplot(status2k.df,cutpoints2k.df,r)
  print(paste0("saving mixplot plot: ",r))
  pdf(paste0(out_dir,'/','GMM_mixture_plot_',r,'.pdf'), width=6.5, height=5)
  print(p)
  dev.off()
}

# save results files
write.csv(status2k.df, file = paste0(out_dir,"/tau_status_2GV.csv"),row.names = F)
write.csv(rates2k.df, file = paste0(out_dir,"/tau_rates_2GV.csv"),row.names = F)
write.csv(cutpoints2k.df, file = paste0(out_dir,"/tau_cutpoints_2GV.csv"),row.names = F)

