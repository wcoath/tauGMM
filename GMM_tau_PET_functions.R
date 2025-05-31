## Functions to run GMM, save pos and rates

# libraries ####
library(tidyverse)
library(mclust)
library(ggpubr)
library(broom)

#extract parameters from GMM
get_params <- function(s) {
  zscore <- 2.326 #99th percentile
  n_k <- s$G
  mod_type=s$modelName
  if (n_k==1) {
    mu1 <- s$mean[1]
    sd1 <- sqrt(s$variance[1])
    mp1 <- s$pro[1]
    #fill in unused
    mu2 <- NA 
    sd2 <- NA
    mp2 <- NA
    mu3 <- NA
    sd3 <- NA
    mp3 <- NA
    cutpoint <- mu1 + zscore * sd1
  }
  else if (n_k==2) {
    mu1 <- s$mean[1]
    sd1 <- sqrt(s$variance[1])
    mp1 <- s$pro[1]
    mu2 <- s$mean[2]
    if (mod_type=="E") {
      #if equal variance, use sd1=sd2
      sd2 <- sd1 }
    else {
      sd2 <- sqrt(s$variance[2])
    }
    mp2 <- s$pro[2]
    mu3 <- NA
    sd3 <- NA
    mp3 <- NA
    if (mu1<mu2) {
      #first gaussian is lower dist
      cutpoint <- mu1 + zscore * sd1
    }
    else {
      #second gaussian is higher dist
      cutpoint <- mu2 + zscore * sd2
    }
  }
  else {
    #three gaussians selected
    mu1 <- s$mean[1]
    sd1 <- sqrt(s$variance[1])
    mp1 <- s$pro[1]
    mu2 <- s$mean[2]
    mp2 <- s$pro[2]
    mu3 <- s$mean[3]
    mp3 <- s$pro[3]
    if (mod_type=="E") {
      
      #if equal variance, use sd1=sd2
      sd2 <- sd1
      sd3 <- sd1 
    }
    else {
      sd2 <- sqrt(s$variance[2])
      sd3 <- sqrt(s$variance[3])
    }
    #select cutpoint based on largest Gaussian. use mixing proportion for this?
    if ( (mp1 >= mp2) & (mp1 >= mp3) ) {
      #first gaussian largest
      cutpoint <- mu1 + zscore * sd1
    }
    else if ( (mp2 >= mp1) & (mp2 >= mp3) )  {
      #second gaussian largest
      cutpoint <- mu2 + zscore * sd2
    }
    else {
      #third gaussian largest
      cutpoint <- mu3 + zscore * sd3
    }
  }
  param_list <- list("n_gaussians" = as.numeric(n_k),"mod_type" = as.character(s$modelName), "cutpoint" = as.numeric(cutpoint),
                     "mu1" = as.numeric(mu1), "mu2" = as.numeric(mu2), "mu3" = as.numeric(mu3), 
                     "sd1" = as.numeric(sd1), "sd2" = as.numeric(sd2), "sd3" = as.numeric(sd3),
                     "mp1" = as.numeric(mp1), "mp2" = as.numeric(mp2), "mp3" = as.numeric(mp3))
  return(param_list)
}

#extract parameters from bootstrapped cutpoints
get_bootcutpoints <- function(bootmod) {
  zscore <- 2.326 #99th percentile
  n_k <- bootmod$G
  mod_type=bootmod$modelName
  if (n_k==1) {
    boot.df <- data.frame("mu1" = bootmod$mean[1:5000,1,1],
                          "sd1" = sqrt(bootmod$variance[1:5000,1,1,1]))
    
    boot.df <- boot.df %>%
      mutate(cutpoint = mu1 + zscore * sd1)
  }
  else if (n_k==2) {
    if (mod_type=="E") {
      boot.df <- data.frame("mu1" = bootmod$mean[1:5000,1,1],
                            "mu2" = bootmod$mean[1:5000,1,2],
                            "sd1" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "sd2" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "mp1" = bootmod$pro[1:5000,1],
                            "mp2" = bootmod$pro[1:5000,2])
    }
    else {
      boot.df <- data.frame("mu1" = bootmod$mean[1:5000,1,1],
                            "mu2" = bootmod$mean[1:5000,1,2],
                            "sd1" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "sd2" = sqrt(bootmod$variance[1:5000,1,1,2]),
                            "mp1" = bootmod$pro[1:5000,1],
                            "mp2" = bootmod$pro[1:5000,2])
    }
    boot.df <- boot.df %>%
      mutate(cutpoint = ifelse(mu1<mu2,
                               mu1 + zscore * sd1,
                               mu2 + zscore * sd2))
  }
  else {
    #three gaussians selected
    if (mod_type=="E") {
      boot.df <- data.frame("mu1" = bootmod$mean[1:5000,1,1],
                            "mu2" = bootmod$mean[1:5000,1,2],
                            "mu3" = bootmod$mean[1:5000,1,3],
                            "sd1" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "sd2" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "sd3" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "mp1" = bootmod$pro[1:5000,1],
                            "mp2" = bootmod$pro[1:5000,2],
                            "mp3" = bootmod$pro[1:5000,3])
    }
    else {
      boot.df <- data.frame("mu1" = bootmod$mean[1:5000,1,1],
                            "mu2" = bootmod$mean[1:5000,1,2],
                            "mu3" = bootmod$mean[1:5000,1,3],
                            "sd1" = sqrt(bootmod$variance[1:5000,1,1,1]),
                            "sd2" = sqrt(bootmod$variance[1:5000,1,1,2]),
                            "sd3" = sqrt(bootmod$variance[1:5000,1,1,3]),
                            "mp1" = bootmod$pro[1:5000,1],
                            "mp2" = bootmod$pro[1:5000,2],
                            "mp3" = bootmod$pro[1:5000,3])
    }
    boot.df <- boot.df %>%
      mutate(cutpoint = ifelse((mp1 >= mp2) & (mp1 >= mp3),
                               mu1 + zscore * sd1,
                               ifelse((mp2 >= mp1) & (mp2 >= mp3),
                                      mu2 + zscore * sd2,
                                      mu3 + zscore * sd3)))
  }
  boot_mean <- mean(boot.df$cutpoint)
  boot_sd <- sd(boot.df$cutpoint)
  boot_low <- as.numeric(quantile(boot.df$cutpoint,0.025))
  boot_high <- as.numeric(quantile(boot.df$cutpoint,0.975))
  boot_list <- list("boot_mean" = boot_mean, "boot_sd" = boot_sd, "boot_low" = boot_low, "boot_high" = boot_high, "boot_cutpoints" = boot.df$cutpoint)
  rm(boot.df)
  return(boot_list)
}

#bootstrap functions, bootstrapping of GMMs to get 5000 cutpoints
#create histogram of bootstrapped cutpoint dist
get_bootplot <- function(bootparams) {
  df <- data.frame("cutpoints" = bootparams$boot_cutpoints)
  p <- ggplot(df, aes(cutpoints)) + 
    geom_histogram(color="black", fill=NA,bins = 100) +
    geom_vline(xintercept=bootparams$boot_mean,color="blue", linetype="dashed", linewidth=1) +
    annotate("rect",xmin=bootparams$boot_low, xmax=bootparams$boot_high, ymin=-Inf, ymax=Inf, alpha=0.2 ,fill="#FF6666") +
    labs(title = 'Bootstrapped GMM Cutpoints (5000 reps)', 
         subtitle = paste0('Mean [95% CI] = ',round(bootparams$boot_mean,digits=3),' [',round(bootparams$boot_low,digits=3),', ',round(bootparams$boot_high,digits = 3),']')) +
    theme_minimal()
  return(p)
}

#remove 5000 cutpoints once summarised in plot and stats
remove_cutpoints <- function(bootparams) {
  boot_list <- list("boot_mean" = bootparams$boot_mean, "boot_sd" = bootparams$boot_sd, "boot_low" = bootparams$boot_low, "boot_high" = bootparams$boot_high)
  return(boot_list)
}


#calculate positivity rates based on cutpoints
get_rates <- function(df) {
  n_pos <- sum(df$status=="Positive")
  n_neg <- sum(df$status=="Negative")
  n_tot <- sum(df$status=="Positive" | df$status=="Negative")
  perc_pos <- (n_pos / n_tot)*100
  perc_neg <- (n_neg / n_tot)*100
  n_pos_low <- sum(df$status_low=="Positive")
  n_neg_low <- sum(df$status_low=="Negative")
  perc_pos_low <- (n_pos_low / n_tot)*100
  perc_neg_low <- (n_neg_low / n_tot)*100
  n_pos_high <- sum(df$status_high=="Positive")
  n_neg_high <- sum(df$status_high=="Negative")
  perc_pos_high <- (n_pos_high / n_tot)*100
  perc_neg_high <- (n_neg_high / n_tot)*100
  
  ratelist <- list("n_pos" = n_pos, "n_neg" = n_neg, "n_tot" = n_tot,
                   "perc_pos" = perc_pos, "perc_neg" = perc_neg, 
                   "n_pos_low" = n_pos_low, "n_neg_low" = n_neg_low,
                   "perc_pos_low" = perc_pos_low, "perc_neg_low" = perc_neg_low,
                   "n_pos_high" = n_pos_high, "n_neg_high" = n_neg_high,
                   "perc_pos_high" = perc_pos_high, "perc_neg_high" = perc_neg_high)
  return(ratelist)
}


get_mixplot <- function(df,cpdf,r) {
  ##generate GMM plot for given method
  # inputs:
  # df = dataframe with suvr values
  # cpdf = dataframe with GMM params
  # r = region label
  
  #need function to plot components:
  plot_mix_comps <- function(x, mu, sigma, lam) {
    # Plot a Mixture Component
    # x = Input data
    # mu = Mean of component
    # sigma = Standard deviation of component
    # lam = Mixture weight of component
    lam * dnorm(x, mu, sigma)
  }
  
  
  #filter dfs for method
  df <- df %>%
    filter(region == r) %>%
    dplyr::select(suvr)
  
  cpdf <- cpdf %>%
    filter(region == r)
  
  #get cutpoint
  cutpoint <- as.numeric(cpdf$cutpoint[1])
  #get number of components
  n_k <- as.character(cpdf$n_gaussians)
  #get model type equal or unequal variance
  mod_type <- as.character(cpdf$mod_type)
  if (mod_type=="E") {
    mod_text<-"equal variances"
  }
  if (mod_type=="V") {
    mod_text<-"unequal variances"
  }
  if (mod_type=="X") {
    mod_text<-"single Gaussian"
  }
  #create plot title adn subtitle
  t <- paste0("Region: ",r,", " ,n_k, " Gaussians with ",mod_text)
  st <- paste0('Cutpoint [95% CI] = ',round(cpdf$cutpoint,digits=3),' [',round(cpdf$boot_low,digits=3),', ',round(cpdf$boot_high,digits = 3),']')
  
  #start plot with first (or only) Gaussian
  p <- df %>%
    ggplot() +
    geom_histogram(aes(suvr, after_stat(density)), binwidth = 0.05,
                   colour = "black", fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(cpdf$mu1, cpdf$sd1, lam = cpdf$mp1),
                  colour = "#2166AC", lwd = 1)
  
  #add extra gaussians
  if (n_k=="2" | n_k=="3") {
    #add extra component if 2 Gaussian
    p <- p +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(cpdf$mu2, cpdf$sd2, lam = cpdf$mp2),
                    colour = "#B2182B", lwd = 1)
  }
  if (n_k=="3") {
    #add extra component if 3 Gaussian
    p <- p +
      stat_function(geom = "line", fun = plot_mix_comps,
                    args = list(cpdf$mu3, cpdf$sd3, lam = cpdf$mp3),
                    colour = "darkolivegreen", lwd = 1)
  }
  
  p <- p +
    labs(x = 'SUVR', y = 'Density', 
         title = t, 
         subtitle = st) +
    geom_vline(xintercept = cutpoint, linetype = 'dashed',col = 'black') +
    annotate("rect", xmin=cpdf$boot_low, xmax=cpdf$boot_high, ymin=-Inf , ymax=Inf, alpha=0.2, fill="red") +
    theme_minimal()
  
  return(p)
}

# Define a function to calculate weighted mean
weighted_mean_region <- function(df, region) {
  df %>%
    dplyr::select(subject,matches(region)) %>%
    rowwise() %>%
    mutate(weighted_mean = weighted.mean(c_across(contains('suvr')), w = c_across(contains('vol')), na.rm = TRUE),
           region = region) %>%
    dplyr::select(subject,region,weighted_mean)
}
