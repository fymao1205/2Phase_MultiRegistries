#ifndef COMMONF_H
#define COMMONF_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

using namespace Rcpp;

double h12_g_A1_X(double t, double X1, String family, NumericVector params0, double beta1, double mu_X2, int logInd);
double survF1_g_A1_X(double t, double X1, String family, NumericVector params0, double beta1, double mu_X2, int logInd);
double h12_g_A1_X_pwc(double t, double X1, NumericVector cuts, NumericVector params0, double beta1, double mu_X2, int logInd);
double survF1_g_A1_X_pwc(double t, double X1, NumericVector cuts, NumericVector params0, double beta1, double mu_X2, int lower, int logInd);
double P12_t_g_A1_X_pwc(double t, double A1, double X1, NumericVector cuts_T2, NumericVector alp_vec, double beta1, double mu_X2_T2, double gamma, NumericVector cuts_TD,NumericVector mu_vec, double eta1, double mu_X2_TD, int logInd);
double P12_t_g_A1_X_pwc_2(double t, double A1, double X1, NumericVector cuts_T2, NumericVector alp_vec, double beta1, double mu_X2_T2, NumericVector cuts_T1D,NumericVector mu1_vec, double eta1, double mu_X2_T1D, NumericVector cuts_T2D,NumericVector mu2_vec, double eta2, double mu_X2_T2D, int logInd);
double P11_t_g_A1_X_pwc(double t, double A1, double X1, NumericVector cuts_T2, NumericVector alp_vec, double beta1, double mu_X2_T2, NumericVector cuts_TD,NumericVector mu_vec, double eta1, double mu_X2_TD, int logInd);
NumericMatrix gaussleg(int n, double x1, double x2);
  
#endif
