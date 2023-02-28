

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

#include "commonf.h"
using namespace Rcpp;

// [[Rcpp::export]]
double logcal_PL_pwc_X1truncexp_num(double Y, double T0, double A1, double T2dgr, double del2,
                                double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                                double mu_X){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  
  double logP_X1_g_X2;
  
  if(X1<1){
    logP_X1_g_X2= R::dexp(X1, 1/mu_X,1);
  }else{
    logP_X1_g_X2= R::pexp(1.0, 1/mu_X, 0, 1);
    }
  
  double lognum = del2*logh12+ logsurvF1+logP_X1_g_X2; //+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD ;
  
  //Rcout << "logP_X1_g_X2 is" << logP_X1_g_X2 << std::endl;
  
  return(lognum);
  
}



// [[Rcpp::export]]
double logcal_PL_pwc_X1truncexp_denom(double Y, double T0, double A1, double T2dgr, double del2,
                                      NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                                  double mu_X){
  
  double res=0.0; 
  
  NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  double denom1=0.0;
  for(int i=0; i<20; i++){
    double tmp = w[i]*survF1_g_A1_X_pwc(T0, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0)*R::dexp(u[i], 1/mu_X,0);
    denom1 = denom1 + tmp;
  }
  double denom11 = survF1_g_A1_X_pwc(T0, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0)*R::pexp(1.0, 1/mu_X,0, 0);
  
  if(Y==1){
    
    //double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (log(denom1+denom11));
    
  }else if(Y==2){
    
    double denom2=1-denom1-denom11;
    
    res = (log(denom2));
    //Rcout << "log denom2 is" << log(denom2) << std::endl;
  } 
  
  return(res);
  
}



// [[Rcpp::export]]
double logobs_PL_pwc_X1truncexp2(double R, double Y, double T0, double A1, double T2dgr, int del2,
                             double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,
                             double mu_X){
  
  double res;
  
  if(R==1){
    
    res = logcal_PL_pwc_X1truncexp_num(Y, T0, A1, T2dgr, del2, 
                                   X1, cuts_T2, params0_T2, beta1, mu_X2_T2, mu_X);
    
  }else{
    
    NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
    NumericVector u = gau_a(_,0);
    NumericVector w = gau_a(_,1);
    
    double res2=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*exp(logcal_PL_pwc_X1truncexp_num(Y, T0, A1, T2dgr, del2, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, mu_X));
      res2 = res2 + tmp; 
    }
    
    double res21 = exp(logcal_PL_pwc_X1truncexp_num(Y, T0, A1, T2dgr, del2, 1.0, cuts_T2, params0_T2, beta1, mu_X2_T2, mu_X));
    
    res= log(res2+res21);
    //if(res2==0){
    //  res = 0;
    //}
  }
  
  double logdenom = logcal_PL_pwc_X1truncexp_denom(Y, T0, A1, T2dgr, del2, cuts_T2, params0_T2, beta1, mu_X2_T2, mu_X);
  
  res = res-logdenom;
  
  
  return(res);
}


// [[Rcpp::export]]
NumericVector obslogPLn_pwc_vec_X1truncexp2(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                        NumericVector T2dgr, IntegerVector del2, 
                                        NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                        NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logobs_PL_pwc_X1truncexp2(R[i], Y[i], T0[i], A1[i], T2dgr[i], del2[i], X1[i],  
                                   cuts_T2, params0_T2, beta1, mu_X2_T2[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logcal_PX1_g_X2_truncexp(double X1, double mu_X){
  
  double logP_X1_g_X2;
  
  if(X1<1){
    logP_X1_g_X2= R::dexp(X1, 1/mu_X,1);
  }else{
    logP_X1_g_X2= R::pexp(1.0, 1/mu_X, 0, 1);
  }
  
  return(logP_X1_g_X2);
  
}

// [[Rcpp::export]]
NumericVector logcaln_PX1_g_X2_vec_X1truncexp2(NumericVector R, NumericVector X1, NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logcal_PX1_g_X2_truncexp(X1[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logcal_PL_pwc_X1truncexp(double Y, double T0, double A1, double T2dgr, double del2,
                                NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  
  double lognum = del2*logh12+ logsurvF1; //+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD ;
  
  //Rcout << "logP_X1_g_X2 is" << logP_X1_g_X2 << std::endl;
  
  //return(lognum);
  
  double res=0.0; 
  
  //NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  //NumericVector u = gau_a(_,0);
  //NumericVector w = gau_a(_,1);
  
  //double denom1=0.0;
  //for(int i=0; i<20; i++){
  //  denom1 =+ w[i]*survF1_g_A1_X_pwc(T0, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0)*R::dexp(u[i], 1/mu_X,0);
  //}
  double denom1 = survF1_g_A1_X_pwc(T0, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0);
  
  if(Y==1){
    
    //double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (lognum-log(denom1));
    
  }else if(Y==2){
    
    double denom2=1-denom1;
    
    res = (lognum-log(denom2));
    //Rcout << "log denom2 is" << log(denom2) << std::endl;
  } 
  
  return(res);
  
  
}


// [[Rcpp::export]]
NumericVector callogPLn_pwc_vec_X1truncexp(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                            NumericVector T2dgr, IntegerVector del2, 
                                            NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2 
                                            ){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logcal_PL_pwc_X1truncexp(Y[i], T0[i], A1[i], T2dgr[i], del2[i],
                                       cuts_T2, params0_T2, beta1, mu_X2_T2[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}



// [[Rcpp::export]]
double logcal_PL_pwc_X1truncexp2(double Y, double T0, double A1, double T2dgr, double del2, double X1,
                                NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2, double mu_X){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  
  double logP_X1_g_X2;
  
  if(X1<1){
    logP_X1_g_X2= R::dexp(X1, 1/mu_X,1);
  }else{
    logP_X1_g_X2= R::pexp(1.0, 1/mu_X, 0, 1);
  }
  
  double lognum = del2*logh12+ logsurvF1+logP_X1_g_X2; //+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD ;
  
  //Rcout << "logP_X1_g_X2 is" << logP_X1_g_X2 << std::endl;
  
  //return(lognum);
  
  double res=0.0; 
  
  //NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  //NumericVector u = gau_a(_,0);
  //NumericVector w = gau_a(_,1);
  
  //double denom1=0.0;
  //for(int i=0; i<20; i++){
  //  denom1 =+ w[i]*survF1_g_A1_X_pwc(T0, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0)*R::dexp(u[i], 1/mu_X,0);
  //}
  double denom1 = survF1_g_A1_X_pwc(T0, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0);
  
  if(Y==1){
    
    //double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (lognum-log(denom1));
    
  }else if(Y==2){
    
    double denom2=1-denom1;
    
    res = (lognum-log(denom2));
    //Rcout << "log denom2 is" << log(denom2) << std::endl;
  } 
  
  return(res);
  
  
}


// [[Rcpp::export]]
NumericVector callogPLn_pwc_vec_X1truncexp2(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                           NumericVector T2dgr, IntegerVector del2, NumericVector X1, 
                                           NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                           NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logcal_PL_pwc_X1truncexp2(Y[i], T0[i], A1[i], T2dgr[i], del2[i], X1[i],
                                      cuts_T2, params0_T2, beta1, mu_X2_T2[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}

// [[Rcpp::export]]
double logcal_PL_pwc_X1truncexp3(double Y, double T0, double A1, double T2dgr, double del2, double X1,
                                 NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2, double mu_X){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  
  double logP_X1_g_X2;
  
  if(X1<1){
    logP_X1_g_X2= R::dexp(X1, 1/mu_X,1);
  }else{
    logP_X1_g_X2= R::pexp(1.0, 1/mu_X, 0, 1);
  }
  
  double lognum = del2*logh12+ logsurvF1+logP_X1_g_X2; //+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD ;
  
  //Rcout << "logP_X1_g_X2 is" << logP_X1_g_X2 << std::endl;
  
  //return(lognum);
  
  double res=0.0; 
  
  NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  double denom1=0.0;
  for(int i=0; i<20; i++){
    double tmp = w[i]*survF1_g_A1_X_pwc(T0, u[i], cuts_T2, params0_T2, 0.0, mu_X2_T2, 0, 0)*R::dexp(u[i], 1/mu_X,0);
    denom1 = denom1 + tmp; 
  }
  double denom11 = survF1_g_A1_X_pwc(T0, 1, cuts_T2, params0_T2, 0.0, mu_X2_T2, 0, 0)*R::pexp(1.0, 1/mu_X,0, 0);
  
  if(Y==1){
    
    //double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (lognum-log(denom1+denom11));
    
  }else if(Y==2){
    
    double denom2=1-denom1-denom11;
    
    res = (lognum-log(denom2));
    //Rcout << "log denom2 is" << log(denom2) << std::endl;
  } 
  
  return(res);
  
  
}


// [[Rcpp::export]]
NumericVector callogPLn_pwc_vec_X1truncexp3(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                            NumericVector T2dgr, IntegerVector del2, NumericVector X1, 
                                            NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                            NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logcal_PL_pwc_X1truncexp3(Y[i], T0[i], A1[i], T2dgr[i], del2[i], X1[i],
                                       cuts_T2, params0_T2, beta1, mu_X2_T2[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}

// [[Rcpp::export]]

double test_gau(double mu_X){
  
  NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  double denom1=0.0;
  for(int i=0; i<20; i++){
    
    double tmp = w[i]*R::dexp(u[i], 1/mu_X,0);
    denom1 = denom1 + tmp;
  }
  
  Rcout << "denom1 is" << denom1 << std::endl;
  
  double denom11 = R::pexp(1.0, 1/mu_X,0, 0);
  
  Rcout << "denom11 is" << denom11 << std::endl;
  
  double res = denom1 + denom11;
  return(res);
  
  }







