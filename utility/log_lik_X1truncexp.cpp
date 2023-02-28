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
double logcal_L_pwc_truncexp2(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                      double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                      NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                      NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D, 
                      double mu_X){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  double logsurvG1_A1 = survF1_g_A1_X_pwc(A1, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 0, 1);
  double logsurvG1_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 0, 1);
  //double logsurvG1_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 0, 1);
  double logsurvG2_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, 0, 1);
  double logsurvG2_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, 0, 1);
  double logeta1D_A2dgr = h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 1);
  double logeta2D_Adgr = h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, 1);
  
  double log_P_obsD = del2*(logsurvG2_Adgr-logsurvG2_A2dgr)+delD*(1-del2)*logeta1D_A2dgr+del2*delD*logeta2D_Adgr;
  
  //double P1_g_X2 = exp(mu_X)/(1+exp(mu_X));
  double logP_X1_g_X2=0.0;
  
  if(X1<1){
    logP_X1_g_X2= R::dexp(X1, 1/mu_X,1);
  }else{
    logP_X1_g_X2= R::pexp(1, 1/mu_X, 0, 1);
  }
  double lognum = del2*logh12+ logsurvF1+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD +logP_X1_g_X2;
  
  //double P11_1, P11_0, P12_1, P12_0, denom1, denom2;
  
  double res=0.0;
  
  NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  
  if(Y==1){
    
    //double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    double denom1=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*R::dexp(u[i], 1/mu_X,0)*P11_t_g_A1_X_pwc(T0, A1, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
      denom1 = denom1+tmp; 
    }
    double denom11 = R::pexp(1, 1/mu_X,0, 0)*P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    
    
    res = (lognum - log(denom1+denom11));
    
  }else if(Y==2){
    
    
    //double P12_1 = P12_t_g_A1_X_pwc_2(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
    //                                  cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    //double P12_0 = P12_t_g_A1_X_pwc_2(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 
    //                                  cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    //double denom2 = P12_1*P1_g_X2 + P12_0*(1-P1_g_X2);  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);      
    
    double denom2=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*R::dexp(u[i], 1/mu_X,0)*P12_t_g_A1_X_pwc_2(T0, A1, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
      denom2 = denom2+tmp; 
    }
    double denom21 = R::pexp(1, 1/mu_X,0, 0)*P12_t_g_A1_X_pwc_2(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    
    res = (lognum - log(denom2+denom21));
    //Rcout << "log denom2 is" << log(denom2) << std::endl;
  } 
  
  //Rcout << "num is" << num << std::endl;
  //Rcout << "lognum is" << lognum << std::endl;
  //Rcout << "P is" << survF1*survG1_A2dgr << std::endl;
  //Rcout << "P11 is" << P11 << std::endl;
  //Rcout << "denom2 is" << denom2 << std::endl;
  
  return(res);
  
}

// [[Rcpp::export]]
double logobs_L_pwc_truncexp2(double R, double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                      double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2, 
                      NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                      NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D, 
                      double mu_X){
  
  double res0, res1, res;
  
  if(R==1){
    
    res = logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                         cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X);
    
  }else{
    //res0 = exp(logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 0.0, cuts_T2, params0_T2, beta1, mu_X2_T2,
    //                          cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    //res1 = exp(logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 1.0, cuts_T2, params0_T2,  beta1, mu_X2_T2,
    //                          cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    
    
    NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
    NumericVector u = gau_a(_,0);
    NumericVector w = gau_a(_,1);
    
    double res0=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*exp(logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2,
                                              cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
      res0 = res0+tmp;
    }
    
    double res1 = exp(logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 1.0, cuts_T2, params0_T2,  beta1, mu_X2_T2,
                                             cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    
    
    double res2 = res0 +res1;
    
    res= log(res2);
    //if(res2==0){
    //  res = 0;
    //}
  }
  
  return(res);
}


// [[Rcpp::export]]
double logcal_L_pwc_truncexp2_num(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                              double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                              NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                              NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D, 
                              double mu_X){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  double logsurvG1_A1 = survF1_g_A1_X_pwc(A1, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 0, 1);
  double logsurvG1_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 0, 1);
  //double logsurvG1_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 0, 1);
  double logsurvG2_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, 0, 1);
  double logsurvG2_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, 0, 1);
  double logeta1D_A2dgr = h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, 1);
  double logeta2D_Adgr = h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, 1);
  
  double log_P_obsD = del2*(logsurvG2_Adgr-logsurvG2_A2dgr)+delD*(1-del2)*logeta1D_A2dgr+del2*delD*logeta2D_Adgr;
  
  //double P1_g_X2 = exp(mu_X)/(1+exp(mu_X));
  double logP_X1_g_X2=0.0;
  
  if(X1<1){
    logP_X1_g_X2= R::dexp(X1, 1/mu_X,1);
  }else{
    logP_X1_g_X2= R::pexp(1, 1/mu_X, 0, 1);
  }
  double lognum = del2*logh12+ logsurvF1+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD +logP_X1_g_X2;
  
  return(lognum);
  
} 



// [[Rcpp::export]]
double logcal_L_pwc_truncexp2_denom(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                              NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                              NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                              NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D, 
                              double mu_X){
  
  double res=0.0;
  
  NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  
  if(Y==1){
    
    //double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    double denom1=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*R::dexp(u[i], 1/mu_X,0)*P11_t_g_A1_X_pwc(T0, A1, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
      denom1 = denom1+tmp; 
    }
    double denom11 = R::pexp(1, 1/mu_X,0, 0)*P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    
    
    res = (log(denom1+denom11));
    
  }else if(Y==2){
    
    
    //double P12_1 = P12_t_g_A1_X_pwc_2(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
    //                                  cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    //double P12_0 = P12_t_g_A1_X_pwc_2(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 
    //                                  cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    //double denom2 = P12_1*P1_g_X2 + P12_0*(1-P1_g_X2);  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);      
    
    double denom2=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*R::dexp(u[i], 1/mu_X,0)*P12_t_g_A1_X_pwc_2(T0, A1, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
      denom2 = denom2+tmp; 
    }
    double denom21 = R::pexp(1, 1/mu_X,0, 0)*P12_t_g_A1_X_pwc_2(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    
    res = (log(denom2+denom21));
    //Rcout << "log denom2 is" << log(denom2) << std::endl;
  } 
  
  //Rcout << "num is" << num << std::endl;
  //Rcout << "lognum is" << lognum << std::endl;
  //Rcout << "P is" << survF1*survG1_A2dgr << std::endl;
  //Rcout << "P11 is" << P11 << std::endl;
  //Rcout << "denom2 is" << denom2 << std::endl;
  
  return(res);
  
}


// [[Rcpp::export]]
double logobs_L_pwc_truncexp3(double R, double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                              double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2, 
                              NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                              NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D, 
                              double mu_X){
  
  double res;
  
  if(R==1){
    
    res = logcal_L_pwc_truncexp2_num(Y, T0, A1, T2dgr, del2, Tdgr, delD, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                                 cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X);
    
  }else{
    //res0 = exp(logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 0.0, cuts_T2, params0_T2, beta1, mu_X2_T2,
    //                          cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    //res1 = exp(logcal_L_pwc_truncexp2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 1.0, cuts_T2, params0_T2,  beta1, mu_X2_T2,
    //                          cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    
    
    NumericMatrix gau_a = gaussleg(20, 0.0, 1.0);
    NumericVector u = gau_a(_,0);
    NumericVector w = gau_a(_,1);
    
    double res0=0.0;
    for(int i=0; i<20; i++){
      double tmp = w[i]*exp(logcal_L_pwc_truncexp2_num(Y, T0, A1, T2dgr, del2, Tdgr, delD, u[i], cuts_T2, params0_T2, beta1, mu_X2_T2,
                                              cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
      res0 = res0 + tmp;
    }
    
    double res1 = exp(logcal_L_pwc_truncexp2_num(Y, T0, A1, T2dgr, del2, Tdgr, delD, 1.0, cuts_T2, params0_T2,  beta1, mu_X2_T2,
                                             cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    
    
    double res2 = res0 +res1;
    
    res= log(res2);
    //if(res2==0){
    //  res = 0;
    //}
  }
  
  double res3 = logcal_L_pwc_truncexp2_denom(Y, T0, A1, T2dgr, del2, Tdgr, delD, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                                           cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X);
  
  return(res-res3);
}



// [[Rcpp::export]]
NumericVector obslogLn_pwc_vec_truncexp2(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                 NumericVector T2dgr, IntegerVector del2, 
                                 NumericVector Tdgr, IntegerVector delD, 
                                 NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                 NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, NumericVector mu_X2_T1D, 
                                 NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, NumericVector mu_X2_T2D, 
                                 NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logobs_L_pwc_truncexp2(R[i], Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i],  
                            cuts_T2, params0_T2, beta1, mu_X2_T2[i], cuts_T1D, params0_T1D, eta1, mu_X2_T1D[i],
                                                                                                           cuts_T2D, params0_T2D, eta2, mu_X2_T2D[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}



// [[Rcpp::export]]
NumericVector obslogLn_pwc_vec_truncexp3(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                         NumericVector T2dgr, IntegerVector del2, 
                                         NumericVector Tdgr, IntegerVector delD, 
                                         NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                         NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, NumericVector mu_X2_T1D, 
                                         NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, NumericVector mu_X2_T2D, 
                                         NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logobs_L_pwc_truncexp3(R[i], Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i],  
                                    cuts_T2, params0_T2, beta1, mu_X2_T2[i], cuts_T1D, params0_T1D, eta1, mu_X2_T1D[i],
                                    cuts_T2D, params0_T2D, eta2, mu_X2_T2D[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}



