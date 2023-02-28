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
double logcal_L_pwc(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                    double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                    double gamma, NumericVector cuts_TD, NumericVector params0_TD, double eta1, double mu_X2_TD, 
                    double mu_X){
  
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  double logsurvG1_A1 = survF1_g_A1_X_pwc(A1, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
  double logsurvG1_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
  double logsurvG1_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
  double logeta1D_A2dgr = h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1);
  double logeta2D_Adgr = h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1)+gamma;
  
  double log_P_obsD = exp(del2*gamma)*(logsurvG1_Adgr-logsurvG1_A2dgr)+delD*(1-del2)*logeta1D_A2dgr+del2*delD*logeta2D_Adgr;
  
  double P1_g_X2 = exp(mu_X)/(1+exp(mu_X));
  double lognum = del2*logh12+ logsurvF1+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD +X1*log(P1_g_X2)+(1-X1)*log(1-P1_g_X2);
  
  //double P11_1, P11_0, P12_1, P12_0, denom1, denom2;
  
  double res=0.0;
  
  if(Y==1){
    
    double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
    double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
    double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (lognum - log(denom1));
    
  }else if(Y==2){
    
    
    double P12_1 = P12_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, gamma, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
    double P12_0 = P12_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, gamma, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
    double denom2 = P12_1*P1_g_X2 + P12_0*(1-P1_g_X2);  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);      
    
    res = (lognum - log(denom2));
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
double logobs_L_pwc(double R, double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                     double X1, 
                     NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2, 
                     double gamma, NumericVector cuts_TD, NumericVector params0_TD, double eta1, double mu_X2_TD, 
                     double mu_X){
  
  double res0, res1, res;
  
  if(R==1){
    
    res = logcal_L_pwc(Y, T0, A1, T2dgr, del2, Tdgr, delD, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                       gamma, cuts_TD, params0_TD, eta1, mu_X2_TD, mu_X);
    
  }else{
    res0 = exp(logcal_L_pwc(Y, T0, A1, T2dgr, del2, Tdgr, delD, 0.0, cuts_T2, params0_T2, beta1, mu_X2_T2,
                            gamma, cuts_TD, params0_TD, eta1, mu_X2_TD, mu_X));
    res1 = exp(logcal_L_pwc(Y, T0, A1, T2dgr, del2, Tdgr, delD, 1.0, cuts_T2, params0_T2,  beta1, mu_X2_T2,
                            gamma, cuts_TD, params0_TD, eta1, mu_X2_TD, mu_X));
    
    double res2 = res0 +res1;
    
    res= log(res2);
    //if(res2==0){
    //  res = 0;
    //}
  }
  
  return(res);
}

// [[Rcpp::export]]
NumericVector obslogLn_pwc_vec(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                                NumericVector T2dgr, IntegerVector del2, 
                                NumericVector Tdgr, IntegerVector delD, 
                                NumericVector X1, //NumericVector X2, 
                                NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                double gamma, NumericVector cuts_TD, NumericVector params0_TD, double eta1, NumericVector mu_X2_TD, 
                                NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logobs_L_pwc(R[i], Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i], //X2[i], 
                          cuts_T2, params0_T2, beta1, mu_X2_T2[i], gamma, cuts_TD, params0_TD, eta1, mu_X2_TD[i], mu_X[i]);
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logcal_L_pwc_2(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
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
  
  double P1_g_X2 = exp(mu_X)/(1+exp(mu_X));
  double lognum = del2*logh12+ logsurvF1+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD +X1*log(P1_g_X2)+(1-X1)*log(1-P1_g_X2);
  
  //double P11_1, P11_0, P12_1, P12_0, denom1, denom2;
  
  double res=0.0;
  
  if(Y==1){
    
    double P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (lognum - log(denom1));
    
  }else if(Y==2){
    
    
    double P12_1 = P12_t_g_A1_X_pwc_2(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                                      cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    double P12_0 = P12_t_g_A1_X_pwc_2(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                                      cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    double denom2 = P12_1*P1_g_X2 + P12_0*(1-P1_g_X2);  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);      
    
    res = (lognum - log(denom2));
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
double logobs_L_pwc_2(double R, double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                    double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2, 
                    NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                    NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D, 
                    double mu_X){
  
  double res0, res1, res;
  
  if(R==1){
    
    res = logcal_L_pwc_2(Y, T0, A1, T2dgr, del2, Tdgr, delD, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                         cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X);
    
  }else{
    res0 = exp(logcal_L_pwc_2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 0.0, cuts_T2, params0_T2, beta1, mu_X2_T2,
                              cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X));
    res1 = exp(logcal_L_pwc_2(Y, T0, A1, T2dgr, del2, Tdgr, delD, 1.0, cuts_T2, params0_T2,  beta1, mu_X2_T2,
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
NumericVector obslogLn_pwc_vec_2(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, 
                               NumericVector T2dgr, IntegerVector del2, 
                               NumericVector Tdgr, IntegerVector delD, 
                               NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                               NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, NumericVector mu_X2_T1D, 
                               NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, NumericVector mu_X2_T2D, 
                               NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    res[i] = logobs_L_pwc_2(R[i], Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i],  
                            cuts_T2, params0_T2, beta1, mu_X2_T2[i], cuts_T1D, params0_T1D, eta1, mu_X2_T1D[i],
                            cuts_T2D, params0_T2D, eta2, mu_X2_T2D[i], mu_X[i]);                                                                                                                                                              
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logcal_PL_pwc(double Y, double T0, double A1, double T2dgr, int del2, double X1,
                     NumericVector cuts, NumericVector params0, double beta1, double mu_X2,
                     double mu_X){
  
  //double logh12, survF1, lognum; 
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 1);
  //double survF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 0, 0);
  
  double P1_g_X2 = exp(mu_X)/(1+exp(mu_X));
  //lognum = del2*logh12 + log(survF1) +X1*log(P1_g_X2)+(1-X1)*log(1-P1_g_X2);
  double lognum = del2*logh12 +survF1_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 0, 1) +X1*log(P1_g_X2)+(1-X1)*log(1-P1_g_X2);
  
  //double survF1_1, survF1_0, denom1, denom2;
  double survF1_1 = survF1_g_A1_X_pwc(T0, 1.0, cuts, params0, beta1, mu_X2, 0, 0);
  double survF1_0 = survF1_g_A1_X_pwc(T0, 0.0, cuts, params0, beta1, mu_X2, 0, 0);
  
  double denom1 = survF1_1*P1_g_X2 + survF1_0*(1-P1_g_X2);
  
  double res=0.0;
  
  if(Y==1){
    
    //Rcout << "lognum is" << lognum << std::endl;
    //Rcout << "denom1 is" << denom1 << std::endl;
    
    //res = (lognum - log(denom1));
    res = lognum - log(denom1);
  
  }else if(Y==2){
    
    double denom2 = 1-denom1;  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);        
    //Rcout << "denom2 is" << denom2 << std::endl;
    res = (lognum - log(denom2));
    
  } 
  
  return(res);
  
}


// [[Rcpp::export]]
double logobs_PL_pwc(double R, double Y, double T0, double A1, double T2dgr, int del2, double X1,
                     NumericVector cuts, NumericVector params0, double beta1, double mu_X2, 
                     double mu_X){
  
  double res;
  
  if(R==1){
    
    res = logcal_PL_pwc(Y, T0, A1, T2dgr, del2, X1, cuts, params0, beta1, mu_X2, mu_X);
    
  }else{
    double res0 = exp(logcal_PL_pwc(Y, T0, A1, T2dgr, del2, 0.0, cuts, params0, beta1, mu_X2, mu_X));
    double res1 = exp(logcal_PL_pwc(Y, T0, A1, T2dgr, del2, 1.0, cuts, params0, beta1, mu_X2, mu_X));
    
    double res2 = res0 +res1;
    
    res= log(res2);
  }
  
  return(res);
}

// [[Rcpp::export]]
NumericVector obslogPLn_pwc_vec(NumericVector R, NumericVector Y, NumericVector T0, NumericVector A1, NumericVector T2dgr, 
                                IntegerVector del2, NumericVector X1, 
                                NumericVector cuts, NumericVector params0, double beta1,
                                NumericVector mu_X2, NumericVector mu_X){
  
  //int n = R.size();
  NumericVector res=no_init(R.size());
  for(int i=0; i<R.size(); i++){
    
    double tmp;
    tmp = logobs_PL_pwc(R[i], Y[i], T0[i], A1[i], T2dgr[i], del2[i],  X1[i], cuts, params0, beta1, mu_X2[i], mu_X[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp);
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logcalL_IPW_pwc(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD, 
                        double X1, //double X2, 
                       NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,
                       double gamma, NumericVector cuts_TD, NumericVector params0_TD, double eta1, double mu_X2_TD){
  
  //double logh12, survF1, etaD_A2dgr, etaD_Adgr, survG1_A2dgr, survG1_Adgr; 
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  double logsurvF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1);
  double logeta1D_A2dgr = h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1);
  double logeta2D_Adgr = h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1)+gamma;
  double logsurvG1_A1 = survF1_g_A1_X_pwc(A1, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
  double logsurvG1_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
  double logsurvG1_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
  
  //double P_obsD = survG1_A2dgr*pow( survG1_Adgr/survG1_A2dgr, exp(del2*gamma))*pow(etaD_A2dgr, delD*(1-del2))*pow(etaD_Adgr, del2*delD)*exp(gamma*del2*delD);
  
  double log_P_obsD = exp(del2*gamma)*(logsurvG1_Adgr-logsurvG1_A2dgr)+delD*(1-del2)*logeta1D_A2dgr+del2*delD*logeta2D_Adgr;
  
  double lognum = del2*logh12 + logsurvF1+(logsurvG1_A2dgr -logsurvG1_A1) + log_P_obsD;
  
  double res=0.0;
  
  if(Y==1){
    
    double P11 = P11_t_g_A1_X_pwc(T0, A1, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_TD, params0_TD, eta1, mu_X2_TD, 0);
    res = (lognum - log(P11));
    
  }else if(Y==2){
    double P12 = P12_t_g_A1_X_pwc(T0, A1, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, gamma, cuts_TD, params0_TD, eta1, mu_X2_TD, 0);
    res = (lognum - log(P12));
    
  } 
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector weighted_obslogLn_IPW_pwc_vec(NumericVector probs, NumericVector Y, NumericVector T0, NumericVector A1, 
                                             NumericVector T2dgr, IntegerVector del2, NumericVector Tdgr, IntegerVector delD, 
                                             NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, 
                                             double beta1, NumericVector mu_X2_T2, double gamma, NumericVector cuts_TD, 
                                             NumericVector params0_TD, double eta1, NumericVector mu_X2_TD){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp= logcalL_IPW_pwc(Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i], //X2[i],
                           cuts_T2, params0_T2, beta1, mu_X2_T2[i],
                           gamma, cuts_TD, params0_TD, eta1, mu_X2_TD[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}

// [[Rcpp::export]]
double logcalL_IPW_pwc_2(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD,
                      double X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,  
                      NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, double mu_X2_T1D, 
                      NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, double mu_X2_T2D){
  
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
  double lognum = del2*logh12+ logsurvF1+(logsurvG1_A2dgr - logsurvG1_A1) + log_P_obsD; //+X1*log(P1_g_X2)+(1-X1)*log(1-P1_g_X2);
  
  //double P11_1, P11_0, P12_1, P12_0, denom1, denom2;
  
  double res=0.0;
  
  if(Y==1){
    
    double P11 = P11_t_g_A1_X_pwc(T0, A1, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_T1D, params0_T1D, eta1, mu_X2_T1D,0);
    //double denom1 = P11_1*P1_g_X2 + P11_0*(1-P1_g_X2);
    //Rcout << "log denom1 is" << log(denom1) << std::endl;
    
    res = (lognum - log(P11));
    
  }else if(Y==2){
    
    
    double P12 = P12_t_g_A1_X_pwc_2(T0, A1, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 
                                      cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    //double P12_0 = P12_t_g_A1_X_pwc_2(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, 
    //                                  cuts_T1D, params0_T1D, eta1, mu_X2_T1D, cuts_T2D, params0_T2D, eta2, mu_X2_T2D,0);
    //double denom2 = P12_1*P1_g_X2 + P12_0*(1-P1_g_X2);  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);      
    
    res = (lognum - log(P12));
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
NumericVector weighted_obslogLn_IPW_pwc_vec2(NumericVector probs, NumericVector Y, NumericVector T0, NumericVector A1, 
                                            NumericVector T2dgr, IntegerVector del2, NumericVector Tdgr, IntegerVector delD, 
                                            NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, double beta1, NumericVector mu_X2_T2, 
                                            NumericVector cuts_T1D, NumericVector params0_T1D, double eta1, NumericVector mu_X2_T1D, 
                                            NumericVector cuts_T2D, NumericVector params0_T2D, double eta2, NumericVector mu_X2_T2D){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp= logcalL_IPW_pwc_2(Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i], 
                                  cuts_T2, params0_T2, beta1, mu_X2_T2[i], cuts_T1D, params0_T1D, eta1, 
                                  mu_X2_T1D[i], cuts_T2D, params0_T2D, eta2, mu_X2_T2D[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logcalPL_T2_pwc(double Y, double T0, double A1, double T2dgr, int del2, double X1, 
                       NumericVector cuts, NumericVector params0, double beta1, double mu_X2){
  
  //double logh12, survF1, lognum; 
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 1);
  double lognum = del2*logh12 + survF1_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 0, 1);
  
  double res=0.0;
  
  if(Y==1){
    //res = (lognum - log(denom));
    
    res = lognum - survF1_g_A1_X_pwc(T0, X1, cuts, params0, beta1, mu_X2, 0, 1);
    
  }else if(Y==2){
    
    double denom = survF1_g_A1_X_pwc(T0, X1, cuts, params0, beta1, mu_X2, 0, 0);     
    res = lognum - log(1-denom);

  } 
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector weighted_obslogPLn_pwc_vec(NumericVector probs, NumericVector Y, NumericVector T0, 
                                         NumericVector A1, NumericVector T2dgr, 
                                         IntegerVector del2, NumericVector X1, 
                                         NumericVector cuts, NumericVector params0, double beta1,
                                         NumericVector mu_X2){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp = logcalPL_T2_pwc(Y[i], T0[i], A1[i], T2dgr[i], del2[i],  X1[i], cuts, params0, beta1, mu_X2[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}

// [[Rcpp::export]]
double logcalPL_T2_wei(double Y, double T0, double A1, double T2dgr, int del2, double X1, 
                       String family, NumericVector params0, double beta1, double mu_X2){
  
  //double logh12, survF1, lognum; 
  double logh12 = h12_g_A1_X(T2dgr, X1, family, params0, beta1, mu_X2, 1);
  //logsurvF1 = survF1_g_A1_X_pwc(T2dgr, A1, X1, cuts, params0, beta1, mu_X2, 0, 1);
  double logsurvF1 = survF1_g_A1_X(T2dgr, X1, family, params0, beta1, mu_X2, 1);
  //double lognum = del2*logh12 +logsurvF1;
  
  double denom= survF1_g_A1_X(T0, X1, family, params0, beta1, mu_X2, 0);         
  
  double res=0.0;
  
  if(Y==1){
    //res = (lognum - log(denom));
    
    res = del2*logh12 +logsurvF1 - log(denom);
    
  }else if(Y==2){
    
    res = del2*logh12 +logsurvF1 - log(1-denom);
    
  } 
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector weighted_obslogPLn_wei_vec(NumericVector probs, NumericVector Y, NumericVector T0, 
                                         NumericVector A1, NumericVector T2dgr, 
                                         IntegerVector del2, NumericVector X1, 
                                         String family, NumericVector params0, double beta1,
                                         NumericVector mu_X2){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp= logcalPL_T2_wei(Y[i], T0[i], A1[i], T2dgr[i], del2[i],  X1[i], family, params0, beta1, mu_X2[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}


// [[Rcpp::export]]
double logscriptPL_T2_wei(double Y, double T0, double A1, double T2dgr, int del2, double X1, 
                          String family, NumericVector params0, double beta1, double mu_X2){
  
  double logh12, survF1, lognum; 
  logh12 = h12_g_A1_X(T2dgr, X1, family, params0, beta1, mu_X2, 1);
  //logsurvF1 = survF1_g_A1_X_pwc(T2dgr, A1, X1, cuts, params0, beta1, mu_X2, 0, 1);
  survF1 = survF1_g_A1_X(T2dgr, X1, family, params0, beta1, mu_X2, 0);
  lognum = del2*logh12 +survF1_g_A1_X(T2dgr, X1, family, params0, beta1, mu_X2, 1);
  
  double res = del2*logh12 +survF1_g_A1_X(T2dgr, X1, family, params0, beta1, mu_X2, 1);
  
  return(res);
  
}


// [[Rcpp::export]]
NumericVector weighted_scriptlogPLn_wei_vec(NumericVector probs, NumericVector Y, NumericVector T0, 
                                         NumericVector A1, NumericVector T2dgr, 
                                         IntegerVector del2, NumericVector X1, 
                                         String family, NumericVector params0, double beta1,
                                         NumericVector mu_X2){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp;
    tmp = logscriptPL_T2_wei(Y[i], T0[i], A1[i], T2dgr[i], del2[i],  X1[i], family, params0, beta1, mu_X2[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}

// [[Rcpp::export]]
double logcalPL_T2_X1_pwc(double Y, double T0, double A1, double T2dgr, int del2, double X1, 
                          NumericVector cuts, NumericVector params0, double beta1, double mu_X2, double p1){
  
  //double logh12, survF1, lognum; 
  double logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 1);
  //logsurvF1 = survF1_g_A1_X_pwc(T2dgr, A1, X1, cuts, params0, beta1, mu_X2, 0, 1);
  //survF1 = survF1_g_A1_X_pwc(T2dgr, A1, X1, cuts, params0, beta1, mu_X2, 0, 0);
  double lognum = del2*logh12 +survF1_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 0, 1)+X1*log(p1)+(1-X1)*log(1-p1);
  
  //double survF1_1, survF1_0, denom1, denom2;
  double survF1_1 = survF1_g_A1_X_pwc(T0, 1.0, cuts, params0, beta1, mu_X2, 0, 0);
  double survF1_0 = survF1_g_A1_X_pwc(T0, 0.0, cuts, params0, beta1, mu_X2, 0, 0);
  
  double denom1 = survF1_1*p1 + survF1_0*(1-p1);
  
  double res=0.0;
  
  if(Y==1){
    
    //Rcout << "lognum is" << lognum << std::endl;
    //Rcout << "denom1 is" << denom1 << std::endl;
    
    //res = (lognum - log(denom1));
    res = lognum - log(denom1);
    
  }else if(Y==2){
    
    double denom2 = 1-denom1;  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);        
    //Rcout << "denom2 is" << denom2 << std::endl;
    res = (lognum - log(denom2));
    
    //if(T0==0){
    //  res=lognum;
    //}
    
  } 
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector weighted_obslogPLn_X1_pwc_vec(NumericVector probs, NumericVector Y, NumericVector T0, 
                                            NumericVector A1, NumericVector T2dgr, 
                                            IntegerVector del2, NumericVector X1, 
                                            NumericVector cuts, NumericVector params0, double beta1,
                                            NumericVector mu_X2, double p1){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp= logcalPL_T2_X1_pwc(Y[i], T0[i], A1[i], T2dgr[i], del2[i],  X1[i], cuts, params0, beta1, mu_X2[i], p1);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}



// [[Rcpp::export]]
double logcalL_IPW_X1_pwc(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD, 
                          double X1, //double X2, 
                          NumericVector cuts_T2, NumericVector params0_T2, double beta1, double mu_X2_T2,
                          double gamma, NumericVector cuts_TD, NumericVector params0_TD, double eta1, double mu_X2_TD, 
                          double p1){
  
  double logh12, survF1, etaD_A2dgr, etaD_Adgr, survG1_A2dgr, survG1_Adgr; 
  logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 1);
  survF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 0);
  etaD_A2dgr = h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0);
  etaD_Adgr = h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0);
  survG1_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 0);///survF1_g_A1_X_pwc(A1, A1, X1, X2, cuts_TD, params0_TD, params_reg_TD);
  survG1_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 0);///survF1_g_A1_X_pwc(A1, A1, X1, X2, cuts_TD, params0_TD, params_reg_TD);
  
  //double P_obsD = survG1_A2dgr*pow( survG1_Adgr/survG1_A2dgr, exp(del2*gamma))*pow(etaD_A2dgr, delD*(1-del2))*pow(etaD_Adgr, del2*delD)*exp(gamma*del2*delD);
  double log_P_obsD = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1)+
    exp(del2*gamma)*(survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1)
                       -survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1))+
                         delD*(1-del2)*h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1)+
                         del2*delD*h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1)+gamma*del2*delD;
  
  double lognum = del2*logh12 + survF1_g_A1_X_pwc(T2dgr, X1, cuts_T2, params0_T2, beta1, mu_X2_T2, 0, 1) + log_P_obsD;
  
  double P11_1, P11_0, P12_1, P12_0, denom1, denom2;
  P11_1 = P11_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
  P11_0 = P11_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
  P12_1 = P12_t_g_A1_X_pwc(T0, A1, 1, cuts_T2, params0_T2, beta1, mu_X2_T2, gamma, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
  P12_0 = P12_t_g_A1_X_pwc(T0, A1, 0, cuts_T2, params0_T2, beta1, mu_X2_T2, gamma, cuts_TD, params0_TD, eta1, mu_X2_TD,0);
  
  denom1 = P11_1*p1 + P11_0*(1-p1);
  denom2 = P12_1*p1 + P12_0*(1-p1);  //(1-survF1_1)*P1_g_X2 + (1-survF1_0)*(1-P1_g_X2);      
  
  double res=0.0;
  
  if(Y==1){
    res = (lognum - log(denom1));
    //if(denom1==0 ){
    //  res=0.0;
    //}
    
  }else if(Y==2){
    
    res = (lognum - log(denom2));
    //if(T0==0){
    //  double survG1_A1 = survF1_g_A1_X_pwc(A1, A1, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 0);///survF1_g_A1_X_pwc(A1, A1, X1, X2, cuts_TD, params0_TD, params_reg_TD);
    //  res=lognum-log(survG1_A1);
    //}
    
  } 
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector weighted_obslogLn_IPW_X1_pwc_vec(NumericVector probs, NumericVector Y, NumericVector T0, NumericVector A1, 
                                               NumericVector T2dgr, IntegerVector del2, NumericVector Tdgr, IntegerVector delD, 
                                               NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, 
                                               double beta1, NumericVector mu_X2_T2, double gamma, NumericVector cuts_TD, 
                                               NumericVector params0_TD, double eta1, NumericVector mu_X2_TD, double p1){
  
  //int n = Y.size();
  NumericVector res(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp;
    tmp = logcalL_IPW_X1_pwc(Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i], 
                             cuts_T2, params0_T2, beta1, mu_X2_T2[i], gamma, cuts_TD, params0_TD, eta1, mu_X2_TD[i], p1);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}



// [[Rcpp::export]]
double logcalPL_T2_TD_pwc(double Y, double T0, double A1, double T2dgr, int del2, double Tdgr, int delD, double X1, 
                          NumericVector cuts, NumericVector params0, double beta1, double mu_X2, 
                          NumericVector cuts_TD, NumericVector params0_TD, double eta1, double mu_X2_TD){
  
  double logh12, survF1, lognum; 
  logh12 = h12_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 1);
  //logsurvF1 = survF1_g_A1_X_pwc(T2dgr, A1, X1, cuts, params0, beta1, mu_X2, 0, 1);
  survF1 = survF1_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 0, 0);
  lognum = del2*logh12 +survF1_g_A1_X_pwc(T2dgr, X1, cuts, params0, beta1, mu_X2, 0, 1);
  
  //double etaD_A2dgr = h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0);
  //double etaD_Adgr = h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0);
  //double survG1_A2dgr = survF1_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 0);///survF1_g_A1_X_pwc(A1, A1, X1, X2, cuts_TD, params0_TD, params_reg_TD);
  //double survG1_Adgr = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 0);///survF1_g_A1_X_pwc(A1, A1, X1, X2, cuts_TD, params0_TD, params_reg_TD);
  
  double log_P_obsD = survF1_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1)+
    delD*(1-del2)*h12_g_A1_X_pwc(A1+T2dgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1)+
    del2*delD*h12_g_A1_X_pwc(A1+Tdgr, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 1);
  
  
  double denom;
  denom = survF1_g_A1_X_pwc(T0, X1, cuts, params0, beta1, mu_X2, 0, 0);         
  
  double res=0.0;
  
  if(Y==1){
    
    res = lognum - survF1_g_A1_X_pwc(T0, X1, cuts, params0, beta1, mu_X2, 0, 1) + log_P_obsD - survF1_g_A1_X_pwc(A1+T0, X1, cuts_TD, params0_TD, eta1, mu_X2_TD, 0, 1);
    
  }else if(Y==2){
    
    res = (lognum - log(1-denom))+ log_P_obsD;
  } 
  
  return(res);
  
}



// [[Rcpp::export]]
NumericVector weighted_obslogPLn_T2_TD_pwc_vec(NumericVector probs, NumericVector Y, NumericVector T0, NumericVector A1, 
                                            NumericVector T2dgr, IntegerVector del2, NumericVector Tdgr, IntegerVector delD, 
                                            NumericVector X1, NumericVector cuts_T2, NumericVector params0_T2, 
                                            double beta1, NumericVector mu_X2_T2, NumericVector cuts_TD, 
                                            NumericVector params0_TD, double eta1, NumericVector mu_X2_TD){
  
  //int n = Y.size();
  NumericVector res=no_init(Y.size());
  for(int i=0; i<Y.size(); i++){
    
    double tmp;
    tmp = logcalPL_T2_TD_pwc(Y[i], T0[i], A1[i], T2dgr[i], del2[i], Tdgr[i], delD[i], X1[i], //X2[i],
                          cuts_T2, params0_T2, beta1, mu_X2_T2[i], cuts_TD, params0_TD, eta1, mu_X2_TD[i]);
    
    //Rcout << "tmp is" << tmp << std::endl;
    
    res[i] = (tmp)/probs[i];
    
  }
  
  return(res);
  
}


