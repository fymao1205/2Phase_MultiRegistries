# 2Phase_MultiRegistries

## Introduction 

This is a repository containing the R code for simulation studies conducted in the project "Two-phase biomarker studies for disease progression with multiple registries". 

## utility

This folder contains all utility functions; 

likelihood function for maximum likelihood estimation and inference purposes are included in *loglik.cpp*, *log_patial_lik_X1truncexp.cpp*, *log_lik_X1truncexp.cpp* and *est_fct.R*; 

Functions for analysis via IPW are included in *ipw_est_fct.R*

*sol_params_fct_general.R*: functions for calculating parameters in different configurations for simulating data sets;

## sim_mle_script.R 

This script is designed to execute an example simulation study within a likelihood-based analysis framework. By modifying the data configurations, it can replicate all the simulation results outlined in sub-Section 4.3 of the manuscript. The script employs a for-loop to conduct repeated simulations. Within each simulation iteration, three key steps are executed: 1) data generation; 2) implementation of the two-phase design; and 3) estimation and inference.

## sim_ipw_script.R 

This script is designed to execute an example simulation study within an IPW framework. By modifying the data configurations, it can replicate all the simulation results outlined in sub-Section 5.3 of the manuscript. The script employs a for-loop to conduct repeated simulations. Within each simulation iteration, three key steps are executed: 1) data generation; 2) implementation of the two-phase design; and 3) estimation and inference.
