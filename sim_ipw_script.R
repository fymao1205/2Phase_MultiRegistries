
# ------------------------------------------------------------------------------------------
# This is a script for one simulation study under IPW analysis
# A for-loop is used for repeated simulations. 
# In each simulation run, 3 steps are included: 
#  1) data generation
#  2) two-phase design implementation
#  3) estimation and inference
# ------------------------------------------------------------------------------------------

library(parallel)
library(survival)
library(MASS)
library(Rcpp)
library(nleqslv)
library(NlcOptim)
source("utility/sol_params_fct_general.R")
source("utility/data_gen.R")
# sourceCpp("utility/commonf.cpp")
# sourceCpp("utility/log_lik.cpp")
# sourceCpp("utility/log_partial_lik_X1truncexp.cpp")
# sourceCpp("utility/log_lik_X1truncexp.cpp")
# source("utility/est_fct.R")
source("utility/ipw_est_fct.R")
source("utility/design_basic.R")
source("utility/design_ipw_suppl.R")


####### Test Data configuration #######

# size of the study population
n=1e6; 
# size of registry 1 in phase I
n1_samp_Y1=1000; 
# size of registry 2 in phase I
n1_samp_Y2=1000; 
# P1: marginal mean of X1; P2: marginal mean of X2; OR: the odds ratio between X1 and X2
P1=0.1; P2=0.5; OR=1
# the birth interval (b0, b1)
b0=0; b1=100
# recruitment date of each disease registry
S0=b1; set_S0=b1; 
# proportion of healthy individuals (in state 0) in the cross-sectional population
pY0=0.5; 
# proportion of individuals in state 1 in the cross-sectional population
pY1=0.01; 
# proportion of individuals in state 2 in the cross-sectional population
pY2=0.01; 
# baseline distribution for all transitions
gamma1=0; gamma2=log(1.5);
family_T1=family_T2=family_TD="weibull"
kai1 =1; kai2 = kaiD = 1.1
# true regression coefficients (X1, X2) for the 0->1 transition
coeffs.beta_T1=c(0, 0.0); 
# true regression coefficients (X1, X2) for transitions to death
coeffs.beta_TD=c(0, 0.0); 
# true regression coefficients (X1, X2, g(T1)) for the disease progression (1->2)
coeffs.beta_T2=c(1, 0.0, 0.0)
# solve for parameters
probvec <- sol.probvec(P1=P1, P2=P2, OR=OR)
pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_gauleg.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                                  family_T1, family_T2, family_TD, 
                                                  kai1, kai2, kaiD, gamma1, gamma2, 
                                                  pY0, pY1, pY2, set_S0)
log.CA <- solve_set_CA_gauleg.f(q_C=0.1, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                family_T1, param_T1=c(exp(pre.lam0.alp0.mu0[1]), kai1), 
                                family_T2, param_T2=c(exp(pre.lam0.alp0.mu0[2]), kai2), 
                                family_TD, param_TD=c(exp(pre.lam0.alp0.mu0[3]), kaiD), gamma1)
lam0=exp(pre.lam0.alp0.mu0[1])
eta0 <- log((1-P2)/probvec[4]-1) 
eta1 <- log(P2/probvec[3]-1)-eta0
eta = c(eta0, eta1)

### user-specified arguments ###

# phase I stratification strategy
phI.design.factor="Y" # e.g., stratigying the phase I data based on registry membership; check design_basics.R to see more options
phI.strat.style = paste0(phI.design.factor, collapse = "_")

# phase II design 
design="bal" # e.g., balanced stratified sampling; check design_basics.R to see more options

# number of simulation runs
nsim=1
# phase II sample size 
n2_samp=900

# user-specified quantiles to determine the cut-points for the PWC hazard functions
# for example, here empirical means of 1->2 and j->3 progression times are used 
q.brks_T2=q.brks_T1D=q.brks_T2D=0.5 
# user-specified quantiles to discretize the continuous stratifying variable for phase I stratification
pcuts=c(0.33,0.66)

# user-specififed to which covariates are included in the primary model of focus, i.e., the 1->2 intensity
# the numbers are corresponding to the column numbers of the covariates matrix
main.cov = c(1,2); 
# user-specififed to which covariates are included in the nuisance covariate model P(X1|X2)
# the numbers are corresponding to the column numbers of the covariates matrix
nui.cov=1

# binary indicator, TRUE: ignore the death transitions in the analysis; FALSE: do not igore 
simplified=FALSE 

iter=1
for(m in 1:(2*nsim)){
  
  if(iter>nsim){break}
  
  set.seed(m)
  
  # generate covariates
  res.X <- gen_x1x2(n, P1, P2, OR) 
  # simulate data from the multistate model 
  dt_0 = gen_illdeath_6state.f(n, X=res.X$data_x1x2, family_01="weibull", params_01=c(lam0, 1), coeffs.beta_01= coeffs.beta_T1,
                               family_12="weibull", params_12=c(exp(pre.lam0.alp0.mu0[2]), 1), coeffs.beta_12= coeffs.beta_T2,
                               family_0D="weibull", params_0D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_0D= coeffs.beta_TD,
                               family_1D="weibull", params_1D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_1D= coeffs.beta_TD,
                               family_2D="weibull", params_2D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_2D= coeffs.beta_TD, 
                               gamma=gamma2)
  # generate data for the study population
  dt_popu = gen_popu.f(n, res.X$data_x1x2, b0, b1, S0, illtime1=dt_0$illtime1, 
                       illtime2=dt_0$illtime2, lifetime=dt_0$deathtime)
  
  # create dataset for disease registry 1
  dt_Y1 = subset(dt_popu, Y1==1 & Y2==0)
  # create dataset for disease registry 2
  dt_Y2 = subset(dt_popu, Y1==1 & Y2==1)
  
  # create the phase I data
  dt_c <- gen_df_obs.f(n1_samp_Y1, n1_samp_Y2, R_Y1=rep(1, n1_samp_Y1), R_Y2=rep(1, n1_samp_Y2), 
                       adcen=exp(log.CA), rcen.param=Inf, dt_Y1, dt_Y2)
  
  # create the covariate matrix 
  Xmat <- as.matrix(cbind(dt_c$X2, log(dt_c$A1_dagger)))
  
  # Phase I stratification 
  phaseI.res <- phase_I_stratify_IDM.f(pcuts, dt_c, phI.strat.style)
  
  # use the user-specified quantiles to determine the cut-points for the PWC models
  brks_T2 <- quantile(subset(dt_c, delta2==1)$T2_dagger, q.brks_T2)
  brks_T1D <- quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T1D)
  brks_T2D <- brks_T1D #quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T2D)
  
  
  # Phase II Selection 
  if(design=="ney"){
    del_beta1 <- try(approx_Delta_beta1_dm.f(gradstep=1e-06, dt_c, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, 
                                             vartheta=c(rep(pre.lam0.alp0.mu0[2], length(brks_T2)+1), 
                                                        c(coeffs.beta_T2), #coeffs.beta_T2[c(2,3,1)],
                                                        rep(pre.lam0.alp0.mu0[1], length(brks_T1D)+1),
                                                        rep(pre.lam0.alp0.mu0[3], length(brks_T2D)+1)
                                             )))
    
    
    s1.res <- Neyman_alloc_influence.f(del_beta1, n2_samp, 
                                       lb.str.n = phaseI.res$strata.n*0.05,
                                       ub.str.n=phaseI.res$strata.n,
                                       phaseI.res)
  }else if(design=="neyA"){
    
    n2_samp_a = n2_samp*prop.a
    s1.res <- Neyman_adpt_2wave_dm.f(gradstep=1e-06, brks_T2, brks_T1D, brks_T2D, n2_samp_a, n2_samp, 
                                     pcuts, phI.design.factor, #phII.design, 
                                     phaseI.res, Xmat, main.cov, simplified)
    
  }else{
    
    s1.res <- subopt.f(n2_samp, design, phaseI.res)
    
  }
  
  # create the dataset after phase II allocation
  dt_s1 <- phaseI.res$dt_ext
  dt_s1[!(dt_s1$id %in% s1.res$s_id), "X1"] <- NA
  dt_s1[!(dt_s1$id %in% s1.res$s_id), "R"] <- 0
  dt_s1$probs <- sapply(dt_s1$group, function(x){
    s1.res$s_prob[x == phaseI.res$ref.strata]
  })
  
  # check how many individuals are selected from each registry
  s.Y = as.data.frame(xtabs(~Y, subset(dt_s1, R==1)))$Freq

  # estimation after a two-phase design
  s1.est <- try(ipw_estL.pwc_dm.f(subset(dt_s1, R==1), Xmat[dt_s1$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
  if(isTRUE(class(s1.est)=="try-error")) { 
    next
  }
  
  # compute standard errors 
  s1.ase <- try(ipw_estL.pwc.ase_dm.f(1e-06, s1.est$par, phaseI.res$ref.strata, subset(dt_s1, R==1), Xmat[dt_s1$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
  if(isTRUE(class(s1.ase)=="try-error")) { 
    next
  }
  
  print(iter)
  iter = iter + 1
  
}

