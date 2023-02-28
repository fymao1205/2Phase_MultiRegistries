
# sourceCpp("~/Dropbox/richard/coding-pj4/illness-death/commonf.cpp")
# sourceCpp("~/Dropbox/richard/coding-pj4/illness-death/log_lik.cpp")

#22/02/17: do not run covariate analysis for mortality intensities
estL.pwc_dm.f <- function(dt, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
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
  }else{
    numpwc_T2D = 1
  }
  
  eta1=eta2=0; 
  mu_X2_T1D = rep(0, dim(Xmat)[1])
  mu_X2_T2D = rep(0, dim(Xmat)[1])
  
  obs.f <- function(vartheta){
    
    params0_T2 = exp(vartheta[1:numpwc_T2])
    beta1 = vartheta[numpwc_T2+1]
    params.reg_T2 = c(vartheta[numpwc_T2+1+(1:length(main.cov))])
    params0_T1D = exp(vartheta[numpwc_T2+length(main.cov)+1+(1:numpwc_T1D)])
    if(simplified){
      params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+2])*params0_T1D
    }else{
      params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
    }
    
    #params.X = c(vartheta[(numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D)+2],0.0)
    #mu_X = as.vector(params.X[1]+outer(Xmat[, nui.cov, drop=F],params.X[-1]))
    #mu_X2_T2 = as.vector(outer(Xmat[, main.cov, drop=F], params.reg_T2))
    #mu_X = as.vector(params.X[1]+(Xmat[, nui.cov, drop=F] %*% params.X[-1]))
    
    mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
    
    if(is.null(nui.cov)){
      mu_X = as.vector(rep(c(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D+2]), dim(Xmat)[1]))
    }else{
      params.X = c(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D+1+(1:(length(nui.cov)+1))])
      mu_X = as.vector(params.X[1]+Xmat[,nui.cov, drop=F] %*% params.X[-1])
    }
    
    res1 = obslogLn_pwc_vec_2(dt$R, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                              dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                              brks_T1D, params0_T1D, eta1, mu_X2_T1D, 
                              brks_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X)
    #res2 = ifelse(is.nan(res1), 0, res1)
    res = -sum(res1)
    res = ifelse(is.finite(res), res, NaN)
    #print(which(is.na(res1)))
    #print(res)
    #print(vartheta)
    #res0 <- scorePL.tr.f(gradstep=1e-06, vartheta, dt)
    
    #res <- colMeans(res0)
    return(res)
  }
  
  #dim <- (numpwc_T2+numpwc_T1D+numpwc_T2D) + length(main.cov)+2 #+length(TD.cov)
  #est.res <- optim(par=rep(-0.01, dim), fn=obs.f, method = c("L-BFGS-B"), hessian = TRUE)
  est.res <- nlm(f=obs.f, p=c(rep(log(0.01), numpwc_T2), rep(0.0, 1+length(main.cov)), 
                              rep(log(0.01), numpwc_T1D), rep(log(0.01), numpwc_T2D), 
                              rep(0, (!is.null(nui.cov))*length(nui.cov)+1)),
                 hessian=T)
  return(est.res)
}

scoreL.pwc_dm.f <- function(gradstep, vartheta, dt, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
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
  }else{
    numpwc_T2D = 1
  } 
  
  eta1=eta2=0; 
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
  #params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
  
  mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
  
  if(is.null(nui.cov)){
    mu_X = as.vector(rep(c(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D+2]), dim(Xmat)[1]))
  }else{
    params.X = c(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D+1+(1:(length(nui.cov)+1))])
    mu_X = as.vector(params.X[1]+Xmat[,nui.cov, drop=F] %*% params.X[-1])
  }

  obsl = obslogLn_pwc_vec_2(dt$R, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                            dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                            brks_T1D, params0_T1D, eta1, mu_X2_T1D, 
                            brks_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X)
  
  dim.param = length(vartheta)
  score.mat = matrix(NA, nrow=(dim(dt)[1]), ncol=dim.param)
  
  for(k in 1:dim.param){
    
    vartheta.p = vartheta
    vartheta.p[k] = vartheta[k] + gradstep
    
    params0_T2 = exp(vartheta.p[1:numpwc_T2])
    beta1 = vartheta.p[numpwc_T2+1]
    params.reg_T2 = c(vartheta.p[numpwc_T2+1+(1:length(main.cov))])
    params0_T1D = exp(vartheta.p[numpwc_T2+length(main.cov)+1+(1:numpwc_T1D)])
    if(simplified){
      params0_T2D = exp(vartheta.p[numpwc_T2+length(main.cov)+numpwc_T1D+2])*params0_T1D
    }else{
      params0_T2D = exp(vartheta.p[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
    }
    #params0_T2D = exp(vartheta.p[numpwc_T2+length(main.cov)+numpwc_T1D+1+(1:numpwc_T2D)])
    
    mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
    
    if(is.null(nui.cov)){
      mu_X = as.vector(rep(c(vartheta.p[numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D+2]), dim(Xmat)[1]))
    }else{
      params.X = c(vartheta.p[numpwc_T2+length(main.cov)+numpwc_T1D+numpwc_T2D+1+(1:(length(nui.cov)+1))])
      mu_X = as.vector(params.X[1]+Xmat[,nui.cov, drop=F] %*% params.X[-1])
    }
    
    obsl.p = obslogLn_pwc_vec_2(dt$R, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                                dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                                brks_T1D, params0_T1D, eta1, mu_X2_T1D, 
                                brks_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X)
    
    score.mat[,k] = (obsl.p - obsl)/gradstep
    
  }
  
  return(score.mat)
}

HessL.pwc_dm.f <- function(gradstep, vartheta, dt, Xmat, main.cov, nui.cov,  brks_T2, brks_T1D, brks_T2D, simplified){
  
  p = length(vartheta)
  
  score_vec_0 <- colSums(scoreL.pwc_dm.f(gradstep, vartheta, dt, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified))
  
  sv_p.list <- lapply(1:p, function(x){
    
    vartheta_p = vartheta 
    vartheta_p[x] = vartheta[x] + gradstep
    
    score_vec_p <- colSums(scoreL.pwc_dm.f(gradstep, vartheta_p, dt, Xmat, main.cov, nui.cov,  brks_T2, brks_T1D, brks_T2D, simplified))
    score_vec_p
  } )
  
  Hess <- matrix(0.0, p, p)
  for(j in 1:p){
    
    Hess[j,] <- (sv_p.list[[j]] - score_vec_0)/gradstep
    
  }
  
  
  return(Hess)
}

sw.aseL.pwc_dm.f <- function(gradstep, vartheta, dt, Xmat, main.cov, nui.cov,  brks_T2, brks_T1D, brks_T2D, simplified){
  
  score.mat = scoreL.pwc_dm.f(gradstep, vartheta, dt, Xmat, main.cov, nui.cov,  brks_T2, brks_T1D, brks_T2D, simplified)
  
  hesmat = HessL.pwc_dm.f(gradstep, vartheta, dt, Xmat, main.cov, nui.cov,  brks_T2, brks_T1D, brks_T2D, simplified)
  
  B = t(score.mat) %*% score.mat
  A = -hesmat
  
  avar = ginv(A)%*% B%*% ginv(A)
  ase = sqrt(diag(avar))
  
  return(ase)
}


