# 2Phase_MultiRegistries

## Introduction 

This is a repository containing the R code for simulation studies conducted in the project "Two-phase biomarker studies for disease progression with multiple registries". 

## utility

This folder contains all utility functions; 

likelihood function for maximum likelihood estimation and inference purposes are included in loglik.cpp, log_patial_lik_X1truncexp.cpp, log_lik_X1truncexp.cpp and est_fct.R; 

Functions for analysis via IPW are included in ipw_est_fct.R

sol_params_fct_general.R: functions for calculating parameters in different configurations for simulating data sets;

## rldt_script.R 

It displays the flow of simulation runs under maximum likelihood and IPW methods. General steps are 1) specifying the parameter settings; 2) solving for unknown parameters; 3) generating a dataset; 4) implement the two-phase design procedure and 5) conduct the estimation adn inference on parameters of interest. 
