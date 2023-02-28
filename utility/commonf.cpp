// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>
  
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double h12_g_A1_X(double t, double X1,  
                  String family, NumericVector params0, double beta1, double mu_X2, int logInd){
  
  double res=0.0;
  //int len_cov = X.size();
  if(family=="weibull"){
    
    double h12_0 = params0[0]*params0[1]*pow(params0[0]*t, params0[1]-1);
    res = h12_0*exp(beta1*X1+mu_X2);
    
    if(logInd) res = log(h12_0) + beta1*X1+mu_X2;
  }
  
  return(res);
}

// [[Rcpp::export]]
double survF1_g_A1_X(double t, double X1, 
                     String family, NumericVector params0, double beta1, double mu_X2, int logInd){
  
  double res=0.0;
  //int len_cov = X.size();
  if(family=="weibull"){
    double cov_reg = exp(beta1*X1+mu_X2);
    double H12_0 = pow(params0[0]*t, params0[1]);
    res = exp(-H12_0*cov_reg);
    if(logInd) 
      res = -H12_0*cov_reg;
  }
  
  return(res);
  
}

//[[Rcpp::export()]]
double hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd)
{
  
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  double y;
  cut.push_front(0);
  cut.push_back(R_PosInf);
  
  if((cut[0] <= x) & (x <= cut[1])){
    y = levels[0];
  }
  if (p > 1.5) {
    for (int i=1; i<p; i++) {
      //y[(cut[i] <= x) & (x <= cut[i + 1])] = levels[i];
      if((cut[i] <= x) & (x <= cut[i+1])){
        y = levels[i];
      }
    }
  }
  if (logInd)
    y = log(y);
  return(y);
}

// [[Rcpp::export]]
double Hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  //NumericVector y(x.size());
  double y;
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  bool which = (cut[0] <= x) & (x <= cut[0 + 1]);
  if (which) {
    //y[which] = x[which];
    y = x*levels[0];
  }
  double wt = levels[0] * cut[1];
  if (p > 1.5) {
    for (int i = 1; i<p; i++) {
      which = (cut[i] <= x) & (x <= cut[i + 1]);
      if (which) {
        //NumericVector xwhich= x[which];
        double tmpx = wt + levels[i] * (x - cut[i]);
        y = tmpx;
      }
      wt = wt + levels[i] * (cut[i + 1] - cut[i]);
    }
  }
  
  //Rcout << "which is " << which << std::endl;
  //Rcout << "wt is " << wt << std::endl;
  //Rcout << "cut is " << cut << std::endl;
  //Rcout << "p is " << p << std::endl;
  //Rcout << "y is " << y << std::endl;
  
  if (logInd)
    y = log(y);
  return(y);
}

//[[Rcpp::export()]]
double ppwc_double(double q, NumericVector cuts, NumericVector levels, int lower, int logInd)
{
  double y;
  if (cuts[0]==0) {
    y = R::pexp(q, 1/levels[0], 0.0, 0.0);
  }else{
    //NumericVector qq(1);
    //qq[0] = q;
    y = Hpwc_double(q, cuts, levels, 0.0);
    if (logInd) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


// [[Rcpp::export]]
double h12_g_A1_X_pwc(double t, double X1, 
                      NumericVector cuts, NumericVector params0, double beta1, double mu_X2, int logInd){
  
  NumericVector levs = params0*exp(beta1*X1+mu_X2);  
  double res = hpwc_double(t, cuts, levs, logInd);
  
  return(res);
}


// [[Rcpp::export]]
double survF1_g_A1_X_pwc(double t, double X1, 
                         NumericVector cuts, NumericVector params0, double beta1, double mu_X2,
                         int lower, int logInd){
  
  NumericVector levs = params0*exp(beta1*X1+mu_X2);  
  double res = ppwc_double(t, cuts, levs, lower, logInd);
  
  return(res);
  
}


// [[Rcpp::export]]
double P11_t_g_A1_X_pwc(double t, double A1, double X1, //X1, double X2,
                        NumericVector cuts_T2, NumericVector alp_vec, double beta1, double mu_X2_T2,  
                        NumericVector cuts_TD, NumericVector mu_vec, double eta1, double mu_X2_TD, int logInd){
  
  //double survF1, survG1;
  //survF1 = survF1_g_A1_X_pwc(t, A1, X1, cuts_T2, alp_vec, beta1, mu_X2_T2, 0, 0);
  //survG1 = survF1_g_A1_X_pwc(A1+t, A1, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 0);//survF1_g_A1_X_pwc(A1, A1, X1, X2, cuts_TD, mu_vec, coeffs_beta_TD);
  
  //Rcout << "res0 is " << res0 << std::endl;
  double res=survF1_g_A1_X_pwc(t, X1, cuts_T2, alp_vec, beta1, mu_X2_T2, 0, 0)*survF1_g_A1_X_pwc(A1+t, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 0)/survF1_g_A1_X_pwc(A1, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 0);
  
  if(logInd)
    res = survF1_g_A1_X_pwc(t, X1, cuts_T2, alp_vec, beta1, mu_X2_T2, 0, 1)+survF1_g_A1_X_pwc(A1+t, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 1)-survF1_g_A1_X_pwc(A1, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 1);
  //double res = lam_X1/(lam_X1+eta0*(1-exp(phi)))*( exp(-eta0*exp(phi)*a) - exp(-(lam_X1+eta0)*a) );
  
  return(res);
  
}

//[[Rcpp::export()]]
NumericMatrix gaussleg(int n, double x1, double x2){
  
  double EPS = 3e-14;
  double m = 0.5*(n+1);
  //double xm, xl;
  double xm = 0.5 * (x2 + x1);
  double xl = 0.5 * (x2 - x1);
  
  //NumericVector x(n), w(n);
  //x = rep(0.0, n);
  //w = rep(0.0, n);
  NumericVector x=no_init(n);
  NumericVector w=no_init(n);
  double z;
  //double tmp;
  for (int i=0; i<m; i++) {
    
    double tmp = M_PI * ( 1.0*(i+1) - 0.25)/(n*1.0 + 0.5);
    z = cos(tmp);
    
    double tol = 9999;
    //double p1, p2, p3;
    double pp;
    //double z1;
    while (tol > EPS) {
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=0; j<n; j++) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2 * j + 1) * z * p2 - j  * p3)/(j+1);
      }
      
      pp = n * (z * p1 - p2)/((z-1)*(1+z));
      double z1 = z;
      z = z1 - p1/pp;
      tol = abs(z - z1);
      //z -= p1/pp;
      //tol= abs(p1/pp);
    }
    
    int s=n -1 - i;
    x[i] = xm - xl * z;
    x[s] = xm + xl * z;
    w[i] = (2 * xl)/((1-z)*(1+z) * pp * pp);
    w[s] = w[i];
  }
  
  NumericMatrix res=no_init_matrix(n, 2);
  res(_,0) = x;
  res(_,1) =w;
  
  //double tmp1 = cos(tmp);
  //Rcout << "tmp is " << tmp << std::endl;
  //Rcout << "tmp1 is " << tmp1 << std::endl;
  return(res);
}


//[[Rcpp::export()]]
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf)
{
  
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  
  y[(cut[0] <= x) & (x < cut[1])] = levels[0];
  if (p > 1.5) {
    for (int i=1; i<p; i++) {
      y[(cut[i] <= x) & (x < cut[i + 1])] = levels[i];
    }
  }
  if (logf)
    y = log(y);
  return(y);
}

//[[Rcpp::export()]]
NumericVector Hpc(NumericVector x,  NumericVector levels, NumericVector cuts, int logf)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  LogicalVector who = (cut[0] <= x) & (x < cut[0 + 1]);
  if (sum(who)) {
    y[who] = x[who];
    y = y*levels[0];
  }
  double su = levels[0] * cut[1];
  if (p > 1.5) {
    for (int i = 1; i<p; i++) {
      who = (cut[i] <= x) & (x < cut[i + 1]);
      if (sum(who)) {
        NumericVector xwho= x[who];
        NumericVector tmpx = su + levels[i] * (xwho - cut[i]);
        y[who] = tmpx;
      }
      su = su + levels[i] * (cut[i + 1] - cut[i]);
    }
  }
  if (logf)
    y = log(y);
  return(y);
}


//[[Rcpp::export()]]
NumericVector vppc(NumericVector q, NumericVector levels,  NumericVector cuts, int lower, int logf)
{
  NumericVector y(1);
  if (cuts[0]==0) {
    y = pexp(q, levels[0], lower, logf);
  }else{
    y = Hpc(q,  levels, cuts, 0.0);
    if (logf) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}

// [[Rcpp::export]]
NumericVector survF1_g_A1_X_pwc_vec(NumericVector t, double X1, 
                                 NumericVector cuts, NumericVector params0, double beta1, double mu_X2,
                                 int lower, int logInd){
  
  NumericVector levs = params0*exp(beta1*X1+mu_X2);  
  NumericVector res = vppc(t, levs, cuts, lower, logInd);
  
  return(res);
  
}

// [[Rcpp::export]]
NumericVector h12_g_A1_X_pwc_vec(NumericVector t, double X1, 
                              NumericVector cuts, NumericVector params0, double beta1, double mu_X2,
                              int logInd){
  
  NumericVector levs = params0*exp(beta1*X1+mu_X2);  
  NumericVector res(t.size());
  res = hpc(t, levs, cuts, logInd) ;
  
  return(res);
  
}


// [[Rcpp::export]]
double P12_t_g_A1_X_pwc(double t, double A1, double X1, 
                     NumericVector cuts_T1, NumericVector alp_vec, 
                     double beta1, double mu_X2_T1, 
                     double gamma, NumericVector cuts_TD, NumericVector mu_vec, 
                     double eta1, double mu_X2_TD, int logInd){
  
  NumericMatrix gau_a = gaussleg(20, 0.0, t);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  NumericVector survF1 = survF1_g_A1_X_pwc_vec(u, X1, cuts_T1, alp_vec, beta1, mu_X2_T1,  0, 0);
  NumericVector logsurvG1 = survF1_g_A1_X_pwc_vec(A1+u, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 1);
  NumericVector h12 = h12_g_A1_X_pwc_vec(u, X1, cuts_T1, alp_vec, beta1, mu_X2_T1, 0);
  NumericVector res0 = survF1*exp(logsurvG1*(1-exp(gamma)))*h12*w;
  double survG1 = survF1_g_A1_X_pwc(A1+t, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 0);
  double res=sum(res0)*pow(survG1, exp(gamma))/survF1_g_A1_X_pwc(A1, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 0);
  
  if(logInd)
    res = exp(gamma)*log(survG1) + log(sum(res0))-survF1_g_A1_X_pwc(A1+t, X1, cuts_TD, mu_vec, eta1, mu_X2_TD, 0, 1);
  
  return(res);
  
}


// [[Rcpp::export]]
double P12_t_g_A1_X_pwc_2(double t, double A1, double X1, 
                        NumericVector cuts_T1, NumericVector alp_vec, double beta1, double mu_X2_T1, 
                        NumericVector cuts_T1D, NumericVector mu1_vec, double eta1, double mu_X2_T1D,
                        NumericVector cuts_T2D, NumericVector mu2_vec, double eta2, double mu_X2_T2D, 
                        int logInd){
  
  NumericMatrix gau_a = gaussleg(20, 0.0, t);
  NumericVector u = gau_a(_,0);
  NumericVector w = gau_a(_,1);
  
  NumericVector survF1 = survF1_g_A1_X_pwc_vec(u, X1, cuts_T1, alp_vec, beta1, mu_X2_T1,  0, 0);
  NumericVector survG1 = survF1_g_A1_X_pwc_vec(A1+u, X1, cuts_T1D, mu1_vec, eta1, mu_X2_T1D, 0, 0);
  NumericVector survG2 = survF1_g_A1_X_pwc_vec(A1+u, X1, cuts_T2D, mu2_vec, eta2, mu_X2_T2D, 0, 0);
  NumericVector h12 = h12_g_A1_X_pwc_vec(u, X1, cuts_T1, alp_vec, beta1, mu_X2_T1, 0);
  NumericVector res0 = (survF1*survG1*h12/survG2)*w;
  double survG2_t = survF1_g_A1_X_pwc(A1+t, X1, cuts_T2D, mu2_vec, eta2, mu_X2_T2D, 0, 0);
  double res=sum(res0)*survG2_t/survF1_g_A1_X_pwc(A1, X1, cuts_T1D, mu1_vec, eta1, mu_X2_T1D, 0, 0);
  
  if(logInd)
    res = log(survG2_t) + log(sum(res0))-survF1_g_A1_X_pwc(A1+t, X1, cuts_T1D, mu1_vec, eta1, mu_X2_T1D, 0, 1);
  
  return(res);
  
}




// [[Rcpp::export]]
NumericVector survF1_g_A1_X2_pwc_vec(NumericVector t, NumericVector A1, 
                                     NumericVector cuts, NumericVector params0, NumericVector mu_X2, int logInd){
  
  //double cov_reg = (params_reg[0]*log(A1) + X1*params_reg[1] + X2*params_reg[2]);
  //double cov_reg = params_reg[0]*A1 + X1*params_reg[1] + X2*params_reg[2];
  int n = t.size();
  NumericVector res(n);
  for(int i=0; i<n; i++){
    //NumericVector levs = params0*exp(mu_X2[i]);  
    //double res0 = ppwc_double(t[i], cuts, levs, 0, 0);
    double res0 = survF1_g_A1_X_pwc(t[i], 1, cuts, params0, 0, mu_X2[i], 0, logInd);
    
    res[i] = res0;
  }
  
  return(res);
  
}


