library(locfit)    #for expit():inverse logit function
library(glue)
library(tidyverse)
library(zoo)  #for na.locf():replacing each NA with the most recent non-NA prior to it
library(patchwork) # for plots patch

#Import the simulation function
source('simulation_function_algorithm1.R')

#generate IPTW results: IPTW estimator's mean, standard deviation, bias and MC standard deviation of bias, save as csv file
#generate corresponding plots for mean & standard deviation of IPTW estimators under different sample size, combine all plots as one plot and save as png file
source('simulation_results_violation_algorithm1.R')




