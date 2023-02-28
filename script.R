library(parallel)
library(survival)
library(MASS)
library(Rcpp)
library(nleqslv)
library(NlcOptim)
source("/work/f5mao/pj4/illness_death/with_X2/DiffMortality/sol_params_fct_general.R")
source("/work/f5mao/pj4/illness_death/data_gen.R")
sourceCpp("/work/f5mao/pj4/illness_death/with_X2/DiffMortality/commonf.cpp")
sourceCpp("/work/f5mao/pj4/illness_death/with_X2/DiffMortality/log_lik.cpp")
source("/work/f5mao/pj4/illness_death/with_X2/DiffMortality/est_fct.R")
source("/work/f5mao/pj4/illness_death/with_X2/design_fct.R")


script_ml_grs_ele_DM.f <- function(nsim, phI.design.factor, phII.design, pcuts, 
                                n, n1_samp_Y1, n1_samp_Y2, n2_samp,
                                P1, P2, OR, b0, b1, S0, 
                                coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                family_T1, family_T2, family_TD, 
                                kai1, kai2, kaiD, gamma1, gamma2, 
                                pY1, pY2, set_S0, q_C, main.cov, nui.cov, q.brks_T2, q.brks_T1D, q.brks_T2D){
  # solve for parameters
  probvec <- sol.probvec(P1=P1, P2=P2, OR=OR)
  
  #pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_exp.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, gamma1, gamma2,
  #                                               pY0, pY1, pY2, set_S0)
  
  #pre.labd_CR <- solve_labdCR_exp.f(q_C, b0, b1, S0, probvec, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
  #                                  lam0=exp(pre.lam0.alp0.mu0[1]), alp0=exp(pre.lam0.alp0.mu0[2]), 
  #                                  mu0=exp(pre.lam0.alp0.mu0[3]), gamma1, pY1)
  
  pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_gauleg.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                 family_T1, family_T2, family_TD, 
                                 kai1, kai2, kaiD, gamma1, gamma2, 
                                 pY0, pY1, pY2, set_S0)
    
  log.CA <- solve_set_CA_gauleg.f(q_C=0.1, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                          family_T1, param_T1=c(exp(pre.lam0.alp0.mu0[1]), kai1), 
                          family_T2, param_T2=c(exp(pre.lam0.alp0.mu0[2]), kai2), 
                          family_TD, param_TD=c(exp(pre.lam0.alp0.mu0[3]), kaiD), gamma1)
  
  lam0=exp(pre.lam0.alp0.mu0[1])
  #labd_CR = exp(pre.labd_CR)
  
  eta0 <- log((1-P2)/probvec[4]-1) #log(p10/(1-P2-p10))
  eta1 <- log(P2/probvec[3]-1)-eta0
  eta = c(eta0, eta1)
  
  #true.par = c(pre.lam0.alp0.mu0[2], coeffs.beta_T2[1:2], gamma2, pre.lam0.alp0.mu0[3], eta)
  
  #print(c(pre.lam0.alp0.mu0, pre.labd_CR))
  
  param.set = c(P1, P2, OR, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, kai1, kai2, kaiD, gamma1, gamma2, pY0, pY1, pY2, q_C)
  name.param.set = c("P1", "P2", "OR", "b1",
                     paste0("phi", 1:length(coeffs.beta_T1)), 
                     paste0("beta", 1:length(coeffs.beta_T2)), 
                     paste0("mu", 1:length(coeffs.beta_TD)),
                     paste0("kai", c("1","2", "D")),
                     "gamma1",  "gamma2","pY0", "pY1", "pY2", "qCR")
  names(param.set) = name.param.set
  
  phI.strat.style = paste0(phI.design.factor, collapse = "_")
  file1.out <- paste(n2_samp, "ML.", phI.strat.style, phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  #file2.out <- paste(n2_samp, "IPW.", phI.strat.style, (length(pcuts)+1), phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  
  if(is.null(q.brks_T2)){
    q.brks_T2 <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T1D)){
    q.brks_T1D <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T2D)){
    q.brks_T2D <- c(0.33, 0.66)
  }
  
  if(phI.strat.style=="Y"){
    ref.group <- 1:2
  }
  
  if(phI.strat.style=="del2"){
    ref.group <- 0:1
  }
  
  if(phI.strat.style=="Y_del2"){
    
    group.mat <- expand.grid(1:2, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    
  }
  
  if(phI.strat.style %in% c("Y_del2_T0", "Y_del2_A1", "Y_del2_T2dgr")){
    
    group.mat <- expand.grid(1:2, 0:1, 1:(length(pcuts)+1))
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    
  }
  
  if(phI.strat.style=="inf_fct"){
    ref.group <- 1:(length(pcuts)+1)
  }
  
  if(phI.strat.style %in% c("del2_T0", "del2_T2dgr", "del2_A1")){
    group.mat <- expand.grid(0:1, 1:(length(pcuts)+1))
    ref.group <- paste0(group.mat[,1], group.mat[,2])
  }
  
  if(is.null(nui.cov)){
    eta.label="eta0"
  }else{
    eta.label <- c("eta0", paste0("eta", nui.cov))
  }
  
  ml.label.tx <- c("nsim", 
                   c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                   paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                   paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                   paste0("eta", c(0, nui.cov))), 
                   paste0("ASE.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                    paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("eta", c(0, nui.cov)))),
                   paste0("score.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                      paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                      paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                                      eta.label)),
                   paste0("g", ref.group), 
                   paste0("sY", 1:2))
  
  
  cat(as.character(ml.label.tx), sep=" ", "\n", append=T, file=file1.out)
  #cat(as.character(ipw.label.tx), sep=" ", "\n", append=T, file=file2.out)
  
  iter=1
  for(m in 1:(2*nsim)){
    
    if(iter>nsim){break}
    
    set.seed(m)
    
    #m_vec = m_mat[m,]
    
    # ------ implement designs to create incomplete dataset
    # generate data
    res.X <- gen_x1x2(n, P1, P2, OR)
    dt_0 = gen_illdeath_6state.f(n, X=res.X$data_x1x2, family_01="weibull", params_01=c(lam0, 1), coeffs.beta_01= coeffs.beta_T1,
                                 family_12="weibull", params_12=c(exp(pre.lam0.alp0.mu0[2]), 1), coeffs.beta_12= coeffs.beta_T2,
                                 family_0D="weibull", params_0D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_0D= coeffs.beta_TD,
                                 family_1D="weibull", params_1D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_1D= coeffs.beta_TD,
                                 family_2D="weibull", params_2D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_2D= coeffs.beta_TD, 
                                 gamma=gamma2)
    
    dt_popu = gen_popu.f(n, res.X$data_x1x2, b0, b1, S0, illtime1=dt_0$illtime1, 
                         illtime2=dt_0$illtime2, lifetime=dt_0$deathtime)
    
    dt_Y1 = subset(dt_popu, Y1==1 & Y2==0)
    dt_Y2 = subset(dt_popu, Y1==1 & Y2==1)
    
    dt_c <- gen_df_obs.f(n1_samp_Y1, n1_samp_Y2, R_Y1=rep(1, n1_samp_Y1), R_Y2=rep(1, n1_samp_Y2), 
                         adcen=exp(log.CA), rcen.param=Inf, dt_Y1, dt_Y2)
    Xmat <- as.matrix(cbind(dt_c$X2, log(dt_c$A1_dagger)))
    phaseI_res <- phase_I_stratify_IDM.f(pcuts, dt_c, phI.strat.style)
    
    brks_T2 <- quantile(subset(dt_c, delta2==1)$T2_dagger, q.brks_T2)
    brks_T1D <- quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T1D)
    brks_T2D <- brks_T1D #quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T2D)
    
    # GRS
    score_res <- try(scoreL_mu_pwc_dm.f(dt_c, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D))
    
    if(isTRUE(class(score_res)=="try-error")) { next }
    
    P_1_g_0 <- mean(dt_c$X1*(1-dt_c$X2))/mean((1-dt_c$X2)) 
    P_1_g_1 <- mean(dt_c$X1*dt_c$X2)/mean(dt_c$X2)
    
    grs.res <- grs_noX1.f(n2_samp, phaseI_res$strata.n, phaseI_res$dt_ext, phaseI_res$ref.strata, 
                          P_1_g_0, P_1_g_1, score_res=score_res)
    dt_grs <- dt_c
    dt_grs[!(dt_grs$id %in% grs.res$s_id), "X1"] <- NA
    dt_grs[!(dt_grs$id %in% grs.res$s_id), "R"] <- 0
    ml.est <- try(estL.pwc_dm.f(dt_grs, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    if(isTRUE(class(ml.est)=="try-error")) { 
      next
    }
    
    ml.ase <- try(sw.aseL.pwc_dm.f(1e-06, ml.est$est, dt_grs, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    if(isTRUE(class(ml.ase)=="try-error")) { 
      next
    }
    
    z <- xtabs(~Y, subset(dt_grs, R==1))
    
    cat(m, ml.est$est, ml.ase, ml.est$gradient, grs.res$s_prob, as.data.frame(z)$Freq,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
  
}

script_ml_rs_ele_DM.f <- function(nsim, phI.design.factor, phII.design, pcuts, 
                                  n, n1_samp_Y1, n1_samp_Y2, n2_samp,
                                  P1, P2, OR, b0, b1, S0, 
                                  coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                  family_T1, family_T2, family_TD, 
                                  kai1, kai2, kaiD, gamma1, gamma2, 
                                  pY1, pY2, set_S0, q_C, main.cov, nui.cov, q.brks_T2, q.brks_T1D, q.brks_T2D){
  # solve for parameters
  probvec <- sol.probvec(P1=P1, P2=P2, OR=OR)
  
  #pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_exp.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, gamma1, gamma2,
  #                                               pY0, pY1, pY2, set_S0)
  
  #pre.labd_CR <- solve_labdCR_exp.f(q_C, b0, b1, S0, probvec, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
  #                                  lam0=exp(pre.lam0.alp0.mu0[1]), alp0=exp(pre.lam0.alp0.mu0[2]), 
  #                                  mu0=exp(pre.lam0.alp0.mu0[3]), gamma1, pY1)
  
  pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_gauleg.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                                    family_T1, family_T2, family_TD, 
                                                    kai1, kai2, kaiD, gamma1, gamma2, 
                                                    pY0, pY1, pY2, set_S0)
  
  log.CA <- solve_set_CA_gauleg.f(q_C=0.1, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                  family_T1, param_T1=c(exp(pre.lam0.alp0.mu0[1]), kai1), 
                                  family_T2, param_T2=c(exp(pre.lam0.alp0.mu0[2]), kai2), 
                                  family_TD, param_TD=c(exp(pre.lam0.alp0.mu0[3]), kaiD), gamma1)
  
  lam0=exp(pre.lam0.alp0.mu0[1])
  #labd_CR = exp(pre.labd_CR)
  
  eta0 <- log((1-P2)/probvec[4]-1) #log(p10/(1-P2-p10))
  eta1 <- log(P2/probvec[3]-1)-eta0
  eta = c(eta0, eta1)
  
  #true.par = c(pre.lam0.alp0.mu0[2], coeffs.beta_T2[1:2], gamma2, pre.lam0.alp0.mu0[3], eta)
  
  #print(c(pre.lam0.alp0.mu0, pre.labd_CR))
  
  param.set = c(P1, P2, OR, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, kai1, kai2, kaiD, gamma1, gamma2, pY0, pY1, pY2, q_C)
  name.param.set = c("P1", "P2", "OR", "b1",
                     paste0("phi", 1:length(coeffs.beta_T1)), 
                     paste0("beta", 1:length(coeffs.beta_T2)), 
                     paste0("mu", 1:length(coeffs.beta_TD)),
                     paste0("kai", c("1","2", "D")),
                     "gamma1",  "gamma2","pY0", "pY1", "pY2", "qCR")
  names(param.set) = name.param.set
  
  phI.strat.style = paste0(phI.design.factor, collapse = "_")
  file1.out <- paste(n2_samp, "ML.", phI.strat.style, phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  #file2.out <- paste(n2_samp, "IPW.", phI.strat.style, (length(pcuts)+1), phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  #if ( file.exists( file2.out ) ) { unlink( file2.out ) }
  
  if(is.null(q.brks_T2)){
    q.brks_T2 <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T1D)){
    q.brks_T1D <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T2D)){
    q.brks_T2D <- c(0.33, 0.66)
  }
  
  if(phI.strat.style=="Y"){
    ref.group <- 1:2
  }
  
  if(phI.strat.style=="del2"){
    ref.group <- 0:1
  }
  
  if(phI.strat.style=="Y_del2"){
    
    group.mat <- expand.grid(1:2, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    
  }
  
  if(phI.strat.style %in% c("Y_del2_T0", "Y_del2_A1", "Y_del2_T2dgr")){
    
    group.mat <- expand.grid(1:2, 0:1, 1:(length(pcuts)+1))
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    
  }
  
  if(phI.strat.style=="inf_fct"){
    ref.group <- 1:(length(pcuts)+1)
  }
  
  if(phI.strat.style %in% c("del2_T0", "del2_T2dgr", "del2_A1")){
    group.mat <- expand.grid(0:1, 1:(length(pcuts)+1))
    ref.group <- paste0(group.mat[,1], group.mat[,2])
  }
  
  if(is.null(nui.cov)){
    eta.label="eta0"
  }else{
    eta.label <- c("eta0", paste0("eta", nui.cov))
  }
  
  ml.label.tx <- c("nsim", 
                   c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                     paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                     paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                     paste0("eta", c(0, nui.cov))), 
                   paste0("ASE.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                    paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("eta", c(0, nui.cov)))),
                   paste0("score.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                      paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                      paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                                      eta.label)),
                   paste0("g", ref.group),
                   paste0("sY", 1:2))
  
  
  cat(as.character(ml.label.tx), sep=" ", "\n", append=T, file=file1.out)
  
  iter=1
  for(m in 1:(2*nsim)){
    
    if(iter>nsim){break}
    
    set.seed(m)
    
    #m_vec = m_mat[m,]
    
    # ------ implement designs to create incomplete dataset
    # generate data
    res.X <- gen_x1x2(n, P1, P2, OR)
    dt_0 = gen_illdeath_6state.f(n, X=res.X$data_x1x2, family_01="weibull", params_01=c(lam0, 1), coeffs.beta_01= coeffs.beta_T1,
                                 family_12="weibull", params_12=c(exp(pre.lam0.alp0.mu0[2]), 1), coeffs.beta_12= coeffs.beta_T2,
                                 family_0D="weibull", params_0D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_0D= coeffs.beta_TD,
                                 family_1D="weibull", params_1D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_1D= coeffs.beta_TD,
                                 family_2D="weibull", params_2D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_2D= coeffs.beta_TD, gamma=gamma2)
    
    dt_popu = gen_popu.f(n, res.X$data_x1x2, b0, b1, S0, illtime1=dt_0$illtime1, 
                         illtime2=dt_0$illtime2, lifetime=dt_0$deathtime)
    
    dt_Y1 = subset(dt_popu, Y1==1 & Y2==0)
    dt_Y2 = subset(dt_popu, Y1==1 & Y2==1)
    
    dt_c <- gen_df_obs.f(n1_samp_Y1, n1_samp_Y2, R_Y1=rep(1, n1_samp_Y1), R_Y2=rep(1, n1_samp_Y2), 
                         adcen=exp(log.CA), rcen.param=Inf, dt_Y1, dt_Y2)
    
    Xmat <- as.matrix(cbind(dt_c$X2, log(dt_c$A1_dagger)))
    
    brks_T2 <- quantile(subset(dt_c, delta2==1)$T2_dagger, q.brks_T2)
    brks_T1D <- quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T1D)
    brks_T2D <- brks_T1D #quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T2D)
    
    
    # Approx. Neyman allocation
    score_res <- try(scoreL_mu_pwc_dm.f(dt_c, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D))
    
    if(isTRUE(class(score_res)=="try-error")) { next }
    
    phaseI_res <- phase_I_stratify_IDM.f(pcuts, dt_c, phI.strat.style)
    
    rs.res <- rs_noX1.f(n2_samp, phaseI_res$strata.n, phaseI_res$dt_ext, phaseI_res$ref.strata, 
                        score_res=score_res)
    dt_rs <- dt_c
    dt_rs[!(dt_rs$id %in% rs.res$s_id), "X1"] <- NA
    dt_rs[!(dt_rs$id %in% rs.res$s_id), "R"] <- 0
    
    ml.est <- try(estL.pwc_dm.f(dt_rs, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    if(isTRUE(class(ml.est)=="try-error")) { 
      next
    }
    
    ml.ase <- try(sw.aseL.pwc_dm.f(1e-06, ml.est$est, dt_rs, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    if(isTRUE(class(ml.ase)=="try-error")) { 
      next
    }
    
    z <- xtabs(~Y, subset(dt_rs, R==1))
    
    cat(m, ml.est$est, ml.ase, ml.est$gradient, rs.res$s_prob, as.data.frame(z)$Freq,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
  
}

script_ml_subopt_ele_DM.f <- function(nsim, phI.design.factor, phII.design, pcuts, 
                                      n, n1_samp_Y1, n1_samp_Y2, n2_samp,
                                      P1, P2, OR, b0, b1, S0, 
                                      coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                      family_T1, family_T2, family_TD, 
                                      kai1, kai2, kaiD, gamma1, gamma2, 
                                      pY1, pY2, set_S0, q_C, main.cov, nui.cov, q.brks_T2, q.brks_T1D, q.brks_T2D){
  
  # solve for parameters
  probvec <- sol.probvec(P1=P1, P2=P2, OR=OR)
  
  #pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_exp.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, gamma1, gamma2,
  #                                               pY0, pY1, pY2, set_S0)
  
  #pre.labd_CR <- solve_labdCR_exp.f(q_C, b0, b1, S0, probvec, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
  #                                  lam0=exp(pre.lam0.alp0.mu0[1]), alp0=exp(pre.lam0.alp0.mu0[2]), 
  #                                  mu0=exp(pre.lam0.alp0.mu0[3]), gamma1, pY1)
  
  pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_gauleg.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                                    family_T1, family_T2, family_TD, 
                                                    kai1, kai2, kaiD, gamma1, gamma2, 
                                                    pY0, pY1, pY2, set_S0)
  
  log.CA <- solve_set_CA_gauleg.f(q_C=0.1, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                  family_T1, param_T1=c(exp(pre.lam0.alp0.mu0[1]), kai1), 
                                  family_T2, param_T2=c(exp(pre.lam0.alp0.mu0[2]), kai2), 
                                  family_TD, param_TD=c(exp(pre.lam0.alp0.mu0[3]), kaiD), gamma1)
  
  lam0=exp(pre.lam0.alp0.mu0[1])
  #labd_CR = exp(pre.labd_CR)
  
  eta0 <- log((1-P2)/probvec[4]-1) #log(p10/(1-P2-p10))
  eta1 <- log(P2/probvec[3]-1)-eta0
  eta = c(eta0, eta1)
  
  #true.par = c(pre.lam0.alp0.mu0[2], coeffs.beta_T2[1:2], gamma2, pre.lam0.alp0.mu0[3], eta)
  
  #print(c(pre.lam0.alp0.mu0, pre.labd_CR))
  
  param.set = c(P1, P2, OR, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, kai1, kai2, kaiD, gamma1, gamma2, pY0, pY1, pY2, q_C)
  name.param.set = c("P1", "P2", "OR", "b1",
                     paste0("phi", 1:length(coeffs.beta_T1)), 
                     paste0("beta", 1:length(coeffs.beta_T2)), 
                     paste0("mu", 1:length(coeffs.beta_TD)),
                     paste0("kai", c("1","2", "D")),
                     "gamma1",  "gamma2","pY0", "pY1", "pY2", "qCR")
  names(param.set) = name.param.set
  
  if(is.null(q.brks_T2)){
    q.brks_T2 <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T1D)){
    q.brks_T1D <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T2D)){
    q.brks_T2D <- c(0.33, 0.66)
  }
  
  phI.strat.style = paste0(phI.design.factor, collapse = "_")
  file1.out <- paste(n2_samp, "ML.", phI.strat.style, phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  #file2.out <- paste(n2_samp, "IPW.", phI.strat.style, (length(pcuts)+1), phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  #if ( file.exists( file2.out ) ) { unlink( file2.out ) }
  
  if(phI.strat.style=="Y"){
    ref.group <- 1:2
  }
  
  if(phI.strat.style=="del2"){
    ref.group <- 0:1
  }
  
  if(phI.strat.style=="Y_del2"){
    
    group.mat <- expand.grid(1:2, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    
  }
  
  if(phI.strat.style %in% c("Y_del2_T0", "Y_del2_A1", "Y_del2_T2dgr")){
    
    group.mat <- expand.grid(1:2, 0:1, 1:(length(pcuts)+1))
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    
  }
  
  if(phI.strat.style=="inf_fct"){
    ref.group <- 1:(length(pcuts)+1)
  }
  
  if(phI.strat.style %in% c("del2_T0", "del2_T2dgr", "del2_A1")){
    group.mat <- expand.grid(0:1, 1:(length(pcuts)+1))
    ref.group <- paste0(group.mat[,1], group.mat[,2])
  }
  
  if(is.null(nui.cov)){
    eta.label="eta0"
  }else{
    eta.label <- c("eta0", paste0("eta", nui.cov))
  }  
  ml.label.tx <- c("nsim", 
                   c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                     paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                     paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                     paste0("eta", c(0, nui.cov))), 
                   paste0("ASE.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                    paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("eta", c(0, nui.cov)))),
                   paste0("score.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                      paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                      paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), 
                                      eta.label)),
                   paste0("g", ref.group),
                   paste0("sY", 1:2))
  
  
  cat(as.character(ml.label.tx), sep=" ", "\n", append=T, file=file1.out)
  #cat(as.character(ipw.label.tx), sep=" ", "\n", append=T, file=file2.out)
  
  iter=1
  for(m in 1:(2*nsim)){
    
    if(iter>nsim){break}
    
    set.seed(m)
    
    #m_vec = m_mat[m,]
    
    # ------ implement designs to create incomplete dataset
    # generate data
    res.X <- gen_x1x2(n, P1, P2, OR)
    dt_0 = gen_illdeath_6state.f(n, X=res.X$data_x1x2, family_01="weibull", params_01=c(lam0, 1), coeffs.beta_01= coeffs.beta_T1,
                                 family_12="weibull", params_12=c(exp(pre.lam0.alp0.mu0[2]), 1), coeffs.beta_12= coeffs.beta_T2,
                                 family_0D="weibull", params_0D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_0D= coeffs.beta_TD,
                                 family_1D="weibull", params_1D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_1D= coeffs.beta_TD,
                                 family_2D="weibull", params_2D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_2D= coeffs.beta_TD, gamma=gamma2)
    
    dt_popu = gen_popu.f(n, res.X$data_x1x2, b0, b1, S0, illtime1=dt_0$illtime1, 
                         illtime2=dt_0$illtime2, lifetime=dt_0$deathtime)
    
    dt_Y1 = subset(dt_popu, Y1==1 & Y2==0)
    dt_Y2 = subset(dt_popu, Y1==1 & Y2==1)
    
    dt_c0 <- gen_df_obs.f(n1_samp_Y1, n1_samp_Y2, R_Y1=rep(1, n1_samp_Y1), R_Y2=rep(1, n1_samp_Y2), 
                          adcen=exp(log.CA), rcen.param=Inf, dt_Y1, dt_Y2)
    
    Xmat <- as.matrix(cbind(dt_c0$X2, log(dt_c0$A1_dagger)))
    
    brks_T2 <- quantile(subset(dt_c0, delta2==1)$T2_dagger, q.brks_T2)
    brks_T1D <- quantile(subset(dt_c0, deltaD==1)$A_dagger, q.brks_T1D)
    brks_T2D <- brks_T1D #quantile(subset(dt_c, deltaD==1)$A_dagger, q.brks_T2D)
    
    
    # Approx. Neyman allocation
    #score_res <- try(scoreL_mu_noX.mortality.f(dt_c0))
    
    #if(isTRUE(class(score_res)=="try-error")) { next }
    #del_mu = c(score_res$M1mu+score_res$aug_S1, score_res$M2mu+score_res$aug_S2)
    
    if(phI.strat.style=="inf_fct"){
      phaseI_res <- phase_I_stratify_Mjmu_IDM.f(pcuts, dt_c0, del_mu)
    }else{
      phaseI_res <- phase_I_stratify_IDM.f(pcuts, dt_c0, phI.strat.style)
    }
    
    # frac0 is specified for # of Y=2
    s1.res <- subopt.f(n2_samp, phII.design, phaseI_res) #ipw_6strata.f(n2_samp, phaseI_res, select_prop=c(rep(2,4),rep(0.05,2)))
    dt_s1 <- phaseI_res$dt_ext
    dt_s1[!(dt_s1$id %in% s1.res$s_id), "X1"] <- NA
    dt_s1[!(dt_s1$id %in% s1.res$s_id), "R"] <- 0
    dt_s1$probs <- sapply(dt_s1$group, function(x){
      s1.res$s_prob[x == phaseI_res$ref.strata]
    })
    
    # ------ estimation and inference
    ml.est <- try(estL.pwc_dm.f(dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    if(isTRUE(class(ml.est)=="try-error")) { 
      next
    }
    
    ml.ase <- try(sw.aseL.pwc_dm.f(1e-06, ml.est$est, dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    if(isTRUE(class(ml.ase)=="try-error")) { 
      next
    }
    
    z <- xtabs(~Y, subset(dt_s1, R==1))
    
    cat(m, ml.est$est, ml.ase, ml.est$gradient, s1.res$s_prob, as.data.frame(z)$Freq,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
  
}


script_ml_srsA_wextMmuB_ele_DM.f <- function(nsim, phI.design.factor, phII.design, pcuts, 
                                          n, n1_samp_Y1, n1_samp_Y2, n2_samp, n2_samp_a,
                                          P1, P2, OR, b0, b1, S0, 
                                          coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                          family_T1, family_T2, family_TD, 
                                          kai1, kai2, kaiD, gamma1, gamma2, 
                                          pY1, pY2, set_S0, q_C, main.cov, nui.cov, 
                                          q.brks_T2, q.brks_T1D, q.brks_T2D){
  # solve for parameters
  probvec <- sol.probvec(P1=P1, P2=P2, OR=OR)
  
  #pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_exp.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, gamma1, gamma2,
  #                                               pY0, pY1, pY2, set_S0)
  
  #pre.labd_CR <- solve_labdCR_exp.f(q_C, b0, b1, S0, probvec, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
  #                                  lam0=exp(pre.lam0.alp0.mu0[1]), alp0=exp(pre.lam0.alp0.mu0[2]), 
  #                                  mu0=exp(pre.lam0.alp0.mu0[3]), gamma1, pY1)
  
  pre.lam0.alp0.mu0 <- solve_lam0_alp0_mu0_gauleg.f(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                                    family_T1, family_T2, family_TD, 
                                                    kai1, kai2, kaiD, gamma1, gamma2, 
                                                    pY0, pY1, pY2, set_S0)
  
  log.CA <- solve_set_CA_gauleg.f(q_C=0.1, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                  family_T1, param_T1=c(exp(pre.lam0.alp0.mu0[1]), kai1), 
                                  family_T2, param_T2=c(exp(pre.lam0.alp0.mu0[2]), kai2), 
                                  family_TD, param_TD=c(exp(pre.lam0.alp0.mu0[3]), kaiD), gamma1)
  
  lam0=exp(pre.lam0.alp0.mu0[1])
  #labd_CR = exp(pre.labd_CR)
  
  eta0 <- log((1-P2)/probvec[4]-1) #log(p10/(1-P2-p10))
  eta1 <- log(P2/probvec[3]-1)-eta0
  eta = c(eta0, eta1)
  
  param.set = c(P1, P2, OR, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, kai1, kai2, kaiD, gamma1, gamma2, pY0, pY1, pY2, q_C)
  name.param.set = c("P1", "P2", "OR", "b1",
                     paste0("phi", 1:length(coeffs.beta_T1)), 
                     paste0("beta", 1:length(coeffs.beta_T2)), 
                     paste0("mu", 1:length(coeffs.beta_TD)),
                     paste0("kai", c("1","2", "D")),
                     "gamma1",  "gamma2","pY0", "pY1", "pY2", "qCR")
  names(param.set) = name.param.set
  
  phI.strat.style = paste0(phI.design.factor, collapse = "_")
  file1.out <- paste(n2_samp, "ML.", phI.strat.style, n2_samp_a, phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  
  if(is.null(q.brks_T2)){
    q.brks_T2 <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T1D)){
    q.brks_T1D <- c(0.33, 0.66)
  }
  
  if(is.null(q.brks_T2D)){
    q.brks_T2D <- c(0.33, 0.66)
  }
  
  if(phI.strat.style=="Y"){
    ref.group <- 1:2
  }
  
  if(phI.strat.style=="del2"){
    ref.group <- 0:1
  }
  
  if(phI.strat.style=="Y_del2"){
    
    group.mat <- expand.grid(1:2, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    
  }
  
  if(phI.strat.style %in% c("Y_del2_T0", "Y_del2_A1", "Y_del2_T2dgr")){
    
    group.mat <- expand.grid(1:2, 0:1, 1:(length(pcuts)+1))
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    
  }
  
  if(phI.strat.style=="inf_fct"){
    ref.group <- 1:(length(pcuts)+1)
  }
  
  if(phI.strat.style %in% c("del2_T0", "del2_T2dgr", "del2_A1")){
    group.mat <- expand.grid(0:1, 1:(length(pcuts)+1))
    ref.group <- paste0(group.mat[,1], group.mat[,2])
  }
  
  if(is.null(nui.cov)){
    eta.label="eta0"
  }else{
    eta.label <- c("eta0", paste0("eta", nui.cov))
  }
  ml.label.tx <- c("nsim", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                             paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                             paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), eta.label),
                   paste0("ASE.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                    paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                    paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), eta.label)), 
                   paste0("score.", c(paste0("logalp.", 1:(length(q.brks_T2)+1)), paste0("beta", c(1, 1+main.cov)), 
                                      paste0("logpsi1.", 1:(length(q.brks_T1D)+1)), 
                                      paste0("logpsi2.", 1:(length(q.brks_T1D)+1)), eta.label)), 
                   paste0("g", ref.group), paste0("s.Y", 1:2))
  
  
  cat(as.character(ml.label.tx), sep=" ", "\n", append=T, file=file1.out)
  
  iter=1
  for(m in 1:(2*nsim)){
    
    if(iter>nsim){break}
    
    set.seed(m)
    
    #m_vec = m_mat[m,]
    
    # ------ implement designs to create incomplete dataset
    # generate data
    res.X <- gen_x1x2(n, P1, P2, OR)
    dt_0 = gen_illdeath_6state.f(n, X=res.X$data_x1x2, family_01="weibull", params_01=c(lam0, 1), coeffs.beta_01= coeffs.beta_T1,
                                 family_12="weibull", params_12=c(exp(pre.lam0.alp0.mu0[2]), 1), coeffs.beta_12= coeffs.beta_T2,
                                 family_0D="weibull", params_0D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_0D= coeffs.beta_TD,
                                 family_1D="weibull", params_1D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_1D= coeffs.beta_TD,
                                 family_2D="weibull", params_2D=c(exp(pre.lam0.alp0.mu0[3]), 1), coeffs.beta_2D= coeffs.beta_TD, gamma=gamma2)
    
    dt_popu = gen_popu.f(n, res.X$data_x1x2, b0, b1, S0, illtime1=dt_0$illtime1, 
                         illtime2=dt_0$illtime2, lifetime=dt_0$deathtime)
    
    dt_Y1 = subset(dt_popu, Y1==1 & Y2==0)
    dt_Y2 = subset(dt_popu, Y1==1 & Y2==1)
    
    dt_c0 <- gen_df_obs.f(n1_samp_Y1, n1_samp_Y2, R_Y1=rep(1, n1_samp_Y1), R_Y2=rep(1, n1_samp_Y2), 
                          adcen=exp(log.CA), rcen.param=Inf, dt_Y1, dt_Y2)
    
    Xmat <- as.matrix(cbind(dt_c0[,c("X2")], log(dt_c0$A1_dagger)))
    
    phaseI_res <- phase_I_stratify_IDM.f(pcuts, dt_c0, phI.strat.style)
    
    brks_T2 <- quantile(subset(dt_c0, delta2==1)$T2_dagger, q.brks_T2)
    brks_T1D <- quantile(subset(dt_c0, deltaD==1)$A_dagger, q.brks_T1D)
    brks_T2D <- brks_T1D #quantile(subset(dt_c0, deltaD==1)$A_dagger, q.brks_T2D)
    
    # Compute Mmu
    Mmu <- try(scoreL_mu_pwc_dm.f(dt_c0, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D))#try(scorePL_mu_pwc.f(phaseI_res$dt_ext, Xmat, main.cov, Inf))
    if(isTRUE(class(Mmu)=="try-error")) { next }
    
    #if(phII.design=="wext_Mmu"){
    
    IIa.res <- subopt.f(n2_samp_a, "srs", phaseI_res) #try(ipw_subopt.f(nsmp_a, design="srs", phaseI_res))
    dt_IIa <- phaseI_res$dt_ext 
    dt_IIa[!(dt_IIa$id %in% IIa.res$s_id), "X1"] <- NA
    dt_IIa[!(dt_IIa$id %in% IIa.res$s_id), "R"] <- 0
    
    phaseI_res.b = phaseI_res
    phaseI_res.b$strata.n <- phaseI_res$strata.n*(1-IIa.res$s_prob)
    phaseI_res.b$Abar_cut <- phaseI_res$Abar_cut
    phaseI_res.b$dt_ext <- subset(phaseI_res$dt_ext, !(id %in% IIa.res$s_id))
    phaseI_res.b$ref.strata <- phaseI_res$ref.strata
    n2_samp_b=n2_samp-n2_samp_a
    
    P_1_g_0 <- mean(subset(dt_IIa,R==1)$X1*(1-subset(dt_IIa,R==1)$X2))/mean((1-subset(dt_IIa,R==1)$X2)) 
    P_1_g_1 <- mean(subset(dt_IIa,R==1)$X1*subset(dt_IIa,R==1)$X2)/mean(subset(dt_IIa,R==1)$X2)
    
    Mmu.b <- Mmu[!(dt_IIa$id %in% IIa.res$s_id)]
    IIb.res <- grs_noX1.f(n2_samp_b, phaseI_res.b$strata.n, phaseI_res.b$dt_ext, phaseI_res.b$ref.strata, 
                          P_1_g_0, P_1_g_1, score_res=Mmu.b)
    
    #}else if(phII.design=="ext_Mmu"){
    
    #  s1.res <- rs_noX1.f(n2_samp, phaseI_res$strata.n, phaseI_res$dt_ext, phaseI_res$ref.strata, 
    #                      score_res=Mmu)
    #}
    
    dt_s1 <- phaseI_res$dt_ext
    dt_s1[!(dt_s1$id %in% c(IIa.res$s_id,IIb.res$s_id)), "X1"] <- NA
    dt_s1[!(dt_s1$id %in% c(IIa.res$s_id,IIb.res$s_id)), "R"] <- 0
    
    ml.est <- try(estL.pwc_dm.f(dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))#try(estPL.pwc.f(dt_s1, Xmat, main.cov, nui.cov, Inf))
    if(isTRUE(class(ml.est)=="try-error")) { 
      next
    }
    ml.ase = try(sw.aseL.pwc_dm.f(1e-06, ml.est$est, dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D)) #sw.asePL.pwc.f(1e-06, ml.est$par, dt_s1, Xmat, main.cov, nui.cov, Inf)
    
    s_m = IIa.res$s_m + IIb.res$s_m
    s_prob = s_m/phaseI_res$strata.n
    
    z <- xtabs(~Y, subset(dt_s1, R==1))
    
    cat(m, ml.est$est, ml.ase, ml.est$gradient, s_prob, as.data.frame(z)$Freq,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
  
}

