

ipw_estL.pwc_dm.f <- function(dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified=TRUE){
  
  if(brks_T2[1]==0){
    numpwc_T2 = 1
  }else{
    numpwc_T2 = 1 + length(brks_T2)
  }
  
  if(brks_T1D[1]==0){
    numpwc_T1D = 1
  }else{
    numpwc_T1D = 1 + length(brks_T1D)
  } 
  
  if(!simplified){
    if(brks_T2D[1]==0){
      numpwc_T2D = 1
    }else{
      numpwc_T2D = 1 + length(brks_T2D)
    } 
  }
  
  mu_X2_T1D = rep(0, dim(Xmat)[1])
  mu_X2_T2D = rep(0, dim(Xmat)[1])
  
  if(simplified){
    
    obs.f <- function(vartheta){
      
      params0_T2 = exp(vartheta[1:numpwc_T2])
      beta1 = vartheta[numpwc_T2+1]
      params.reg_T2 = c(vartheta[numpwc_T2+1+(1:length(main.cov))])
      params0_T1D = exp(vartheta[numpwc_T2+length(main.cov)+1+(1:numpwc_T1D)])
      params0_T2D = params0_T1D*exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+2]) #exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
      
      mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
      
      res1 = weighted_obslogLn_IPW_pwc_vec2(dt$probs, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                                            dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                                            brks_T1D, params0_T1D, 0, mu_X2_T1D, 
                                            brks_T2D, params0_T2D, 0, mu_X2_T2D)
      res1 <- ifelse(is.finite(res1), res1, 0)
      res = -sum(res1)
      
      #print(res)
      #print(vartheta)
      #obj.res <- ipw_weighted_Pscore_noX2.pwc.f(1e-06, dt, vartheta, brks)
      #res <- colSums(obj.res) 
      
      return(res)
    }
    
    dim <- numpwc_T2+numpwc_T1D+2+length(main.cov)
    #est.res <- nlm(obs.f, p=rep(0.05, dim), gradtol=1e-07, steptol=1e-07, hessian=TRUE)
    est.res <- optim(par=rep(-0.01, dim), 
                     fn=obs.f, method = c("L-BFGS-B"), hessian = TRUE)
    #est.res <- nlm(f=obs.f, p=c(rep(log(0.01), numpwc_T2), rep(0.0, 1+length(main.cov)), 
    #                            rep(log(0.01), numpwc_T1D), rep(log(0.01), numpwc_T2D)),
    #               hessian=T)
    
  }else{
    
    obs.f <- function(vartheta){
      
      params0_T2 = exp(vartheta[1:numpwc_T2])
      beta1 = vartheta[numpwc_T2+1]
      params.reg_T2 = c(vartheta[numpwc_T2+1+(1:length(main.cov))])
      params0_T1D = exp(vartheta[numpwc_T2+length(main.cov)+1+(1:numpwc_T1D)])
      params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
      
      mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
      
      res1 = weighted_obslogLn_IPW_pwc_vec2(dt$probs, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                                            dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                                            brks_T1D, params0_T1D, 0, mu_X2_T1D, 
                                            brks_T2D, params0_T2D, 0, mu_X2_T2D)
      res1 <- ifelse(is.finite(res1), res1, 0)
      res = -sum(res1)
      
      #print(res)
      #print(vartheta)
      #obj.res <- ipw_weighted_Pscore_noX2.pwc.f(1e-06, dt, vartheta, brks)
      #res <- colSums(obj.res) 
      
      return(res)
    }
    
    dim <- numpwc_T2+numpwc_T1D+numpwc_T2D+1+length(main.cov)
    #est.res <- nlm(obs.f, p=rep(0.05, dim), gradtol=1e-07, steptol=1e-07, hessian=TRUE)
    est.res <- optim(par=rep(-0.01, dim), 
                     fn=obs.f, method = c("L-BFGS-B"), hessian = TRUE)
    #est.res <- nlm(f=obs.f, p=c(rep(log(0.01), numpwc_T2), rep(0.0, 1+length(main.cov)), 
    #                            rep(log(0.01), numpwc_T1D), rep(log(0.01), numpwc_T2D)),
    #               hessian=T)
    
  }
  
  
  return(est.res)
}

ipw_weighted_score.pwc_dm.f <- function(gradstep, vartheta, dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified=TRUE){
  
  if(brks_T2[1]==0){
    numpwc_T2 = 1
  }else{
    numpwc_T2 = 1 + length(brks_T2)
  }
  
  if(brks_T1D[1]==0){
    numpwc_T1D = 1
  }else{
    numpwc_T1D = 1 + length(brks_T1D)
  } 
  
  if(!simplified){
    if(brks_T2D[1]==0){
      numpwc_T2D = 1
    }else{
      numpwc_T2D = 1 + length(brks_T2D)
    } 
  }
  
  #eta1=eta2=0; 
  mu_X2_T1D = rep(0, dim(Xmat)[1])
  mu_X2_T2D = rep(0, dim(Xmat)[1])
  
  params0_T2 = exp(vartheta[1:numpwc_T2])
  beta1 = vartheta[numpwc_T2+1]
  params.reg_T2 = c(vartheta[numpwc_T2+1+(1:length(main.cov))])
  params0_T1D = exp(vartheta[numpwc_T2+length(main.cov)+1+(1:numpwc_T1D)])
  
  if(simplified){
    params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+2])*params0_T1D
  }else{
    params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
  }
  
  mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
  
  logL_0 = weighted_obslogLn_IPW_pwc_vec2(dt$probs, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                                        dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                                        brks_T1D, params0_T1D, 0, mu_X2_T1D, 
                                        brks_T2D, params0_T2D, 0, mu_X2_T2D)
  
  dim.param = length(vartheta)
  score.mat = matrix(NA, nrow=(dim(dt)[1]), ncol=dim.param)
  
  score.list <- lapply(1:dim.param, function(x){
    
    vartheta_p = vartheta 
    vartheta_p[x] = vartheta[x] + gradstep
    
    params0_T2 = exp(vartheta_p[1:numpwc_T2])
    beta1 = vartheta_p[numpwc_T2+1]
    params.reg_T2 = c(vartheta_p[numpwc_T2+1+(1:length(main.cov))])
    params0_T1D = exp(vartheta_p[numpwc_T2+length(main.cov)+1+(1:numpwc_T1D)])
    #params0_T2D = exp(vartheta_p[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
    
    if(simplified){
      params0_T2D = exp(vartheta_p[numpwc_T2+length(main.cov)+numpwc_T1D+2])*params0_T1D
    }else{
      params0_T2D = exp(vartheta_p[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
    }
    
    mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
    
    logL_p <- weighted_obslogLn_IPW_pwc_vec2(dt$probs, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                                             dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                                             brks_T1D, params0_T1D, 0, mu_X2_T1D, 
                                             brks_T2D, params0_T2D, 0, mu_X2_T2D)
    
    tmp <- (logL_p - logL_0)/gradstep
    tmp
  } )
  
  score_mat <- do.call("cbind", score.list)
  
  return(score_mat)
}

ipw_weighted_hess.pwc_dm.f <- function(gradstep, vartheta, dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
  p = length(vartheta)
  
  score_vec_0 <- colSums(ipw_weighted_score.pwc_dm.f(gradstep, vartheta, dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
  
  sv_p.list <- lapply(1:p, function(x){
    
    vartheta_p = vartheta 
    vartheta_p[x] = vartheta[x] + gradstep
    
    score_vec_p <- colSums(ipw_weighted_score.pwc_dm.f(gradstep, vartheta_p, dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
    score_vec_p
  } )
  
  Hess <- matrix(0.0, p, p)
  for(j in 1:p){
    
    Hess[j,] <- (sv_p.list[[j]] - score_vec_0)/gradstep
    
  }
  
  
  return(Hess)
}

ipw_estL.pwc.ase_dm.f <- function(gradstep, vartheta, ref.group, dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
  # score_vec
  #score_mat <- ipw_weighted_Pscore.pwc.f(gradstep, subset(dt_use, R==1), est.res$x, brks) # est.res$par
  score_mat <- ipw_weighted_score.pwc_dm.f(gradstep, vartheta, subset(dt, R==1), Xmat[dt$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, simplified)
  #score_mat <- ipw_weighted_Pscore.tr.f(gradstep, subset(dt_use, R==1), vartheta)
  B11 = t(score_mat) %*% score_mat
  
  # Hess
  #A11 = -est.res$hessian
  A11 = -ipw_weighted_hess.pwc_dm.f(gradstep, vartheta, subset(dt, R==1), Xmat[dt$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, simplified)
  
  # A12
  score_probs_mat <- sapply(1:(dim(subset(dt, R==1))[1]), function(i){
    ind = (subset(dt, R==1)$group[i] == ref.group)
    res = numeric(length(ref.group))
    res[ind] = 1/subset(dt, R==1)$probs[i]
    res
  })
  score_probs_mat = t(score_probs_mat)
  A12 = t(score_mat) %*% score_probs_mat
  
  # A22
  hess_probs <- sapply(1:(dim(dt)[1]), function(i){
    ind = (dt$group[i] == ref.group)
    res = numeric(length(ref.group))
    num = dt$R[i] - dt$probs[i]
    denom = dt$probs[i]*(1-dt$probs[i])
    res[ind] = ifelse(denom==0, 0, (num/denom)^2)
    diag(res)
  })
  A22 <-matrix(rowSums(hess_probs), nrow=length(ref.group), byrow=0)
  pos.nonzero <- which(diag(A22)!=0)
  if(any(diag(A22)!=0)){
    A22 <- A22[pos.nonzero, pos.nonzero]
    A12 <- A12[,pos.nonzero]
  }
  
  avar <- ginv(A11) %*%(B11 - A12 %*% ginv(A22) %*% t(A12)) %*% t(ginv(A11))
  ase = sqrt(diag(avar))
  
  return(list(se=ase, score=colMeans(score_mat)))
}

