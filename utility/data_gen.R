# data generation

### generate (X1,X2)
# P(X1)=0.5, P(X2)=0.5, or 0.2, OR=2
# With margins P(X1=1)=P1 and P(X2=1)=P2, and Odds ratio OR, generate (x1,x2) 

logit.f <- function(x){
  res <- log(x/(1-x))
  return(res)
}

expit.f <- function(x){
  
  exp(x)/(1+exp(x))
}

P_X1_g_X2.f <- function(X, param.cov){
  
  if(is.null(dim(X))){
    
    X1 <- X[1]
    X2 <- X[2]
    
  }else{
    
    X1 <- X[,1]
    X2 <- X[,2]
    
  }
  
  res1 <- expit.f(param.cov[1]+param.cov[2]*X2)
  res <- res1^X1*(1-res1)^(1-X1)
  return(res)
}

gen_x1x2 <- function(n, P1, P2, OR)
{
  
  probvec <- sol.probvec(P1,P2,OR)
  
  #p00 = probvec[1]
  p10 = probvec[2]
  eta0 = log(p10/(1-P2-p10))
  eta1 = log(OR)
  x2 <- rbinom(n, 1, P2)
  z <- eta0 + eta1*x2 
  pr = expit.f(z)
  x1 <- rbinom(n, 1, pr)
  
  data_x1x2 <- cbind(x1,x2)
  true.eta <- c(eta0, eta1)
  
  return(list(data_x1x2=data_x1x2,
              true.eta=true.eta)) 
  
}

### generate transition times given X

T.f <- function(n, X, coeffs.beta, params, family="weibull"){
  
  u_n <- runif(n)
  
  if(family=="weibull" | is.null(family)){
    scale_pm <- params[1]
    shape_pm <- params[2]
    # t_n <- qweibull(1- (1-u_n)^(exp(-X%*%coeffs.beta)),
    #                 shape=shape_pm, scale=1/scale_pm)
    res1 <- -log(1-u_n)*exp(-X%*%coeffs.beta)
    t_n <- res1^(1/shape_pm)/scale_pm
  }
  
  if(family=="pwc"){
    dim_pm <- length(params)
    pwc_vec <- params[1:((dim_pm+1)/2)]
    brks_vec <- params[-(1:((dim_pm+1)/2))]
    t_n <- qpch(1- u_n^(exp(-X%*%coeffs.beta)), cuts=brks_vec, levels=pwc_vec)
  }
  
  return(t_n)
  
}

T_1D.f <- function(n, X, coeffs.beta, gamma, params, family="weibull"){
  
  u_n <- runif(n)
  
  if(family=="weibull" | is.null(family)){
    scale_pm <- params[1]
    shape_pm <- params[2]
    # t_n <- qweibull(1- (1-u_n)^(exp(-X%*%coeffs.beta)),
    #                 shape=shape_pm, scale=1/scale_pm)
    res1 <- -log(1-u_n)*exp(-X%*%coeffs.beta-gamma)
    t_n <- res1^(1/shape_pm)/scale_pm
  }
  
  if(family=="pwc"){
    dim_pm <- length(params)
    pwc_vec <- params[1:((dim_pm+1)/2)]
    brks_vec <- params[-(1:((dim_pm+1)/2))]
    t_n <- qpch(1- u_n^(exp(-X%*%coeffs.beta-gamma)), cuts=brks_vec, levels=pwc_vec)
  }
  
  return(t_n)
  
}

T2_g_T1.f <- function(n, T1, X, coeffs.beta, params,  family="weibull"){
  
  u_n <- runif(n)
  
  if(family=="weibull" | is.null(family)){
    scale_pm <- params[1]
    shape_pm <- params[2]
    # t_n <- qweibull(1- (1-u_n)^(exp(-X%*%coeffs.beta)),
    #                 shape=shape_pm, scale=1/scale_pm)
    res1 <- -log(1-u_n)*exp(-X%*%coeffs.beta[-1] - log(T1)*coeffs.beta[1])
    t_n <- res1^(1/shape_pm)/scale_pm
  }
  
  if(family=="pwc"){
    dim_pm <- length(params)
    pwc_vec <- params[1:((dim_pm+1)/2)]
    brks_vec <- params[-(1:((dim_pm+1)/2))]
    t_n <- qpch(1- u_n^(exp(-X%*%coeffs.beta[-1] - log(T1)*coeffs.beta[1])), cuts=brks_vec, levels=pwc_vec)
  }
  
  return(t_n)
  
}

#n=1e6; 
#X <- gen_x1x2(n, P1=0.5, P2=0.5, OR=2)$data_x1x2
#params_01 <- c(0.01,1); params_12 = c(0.02, 1)
#params_0D <- c(0.006, 1); params_1D <- c(0.01, 1); params_2D <- c(0.01, 1)
#family_01 <- family_12 <- family_2D <- family_1D <- family_0D <- "weibull"
#coeffs.beta_01 <- c(log(1), 0.25); coeffs.beta_12 <- c(log(1.1), log(1), 0.25) 
#coeffs.beta_0D <- c(0, log(1.1)); coeffs.beta_1D <- coeffs.beta_2D <- c(0, log(1.1))
#gamma=0


gen_illdeath_6state.f <- function(n, X,   
                           family_01="weibull", params_01, coeffs.beta_01,
                           family_12="weibull", params_12, coeffs.beta_12,
                           family_0D="weibull", params_0D, coeffs.beta_0D,
                           family_1D="weibull", params_1D, coeffs.beta_1D,
                           family_2D="weibull", params_2D, coeffs.beta_2D, gamma){
  
  
  ### generate T_01
  t_01_n <- T.f(n,X, coeffs.beta_01, params_01, family=family_01)
  
  ### generate T_1D
  t_1D_n <- T_1D.f(n,X, coeffs.beta_1D, gamma=0, params_1D, family=family_1D)
  
  ### generate T_2D
  t_2D_n <- T_1D.f(n,X, coeffs.beta_2D, gamma, params_2D, family=family_2D)
  
  ### generate T_0D
  t_0D_n <- T.f(n,X, coeffs.beta_0D, params_0D, family=family_0D)
  
  ### generate T_12
  delta_01 <- 1*(t_01_n <= t_0D_n)
  delta_0D <- 1-delta_01
  
  t_12_n <- T2_g_T1.f(n, t_01_n, X, coeffs.beta_12, params_12, family_12)
  
  delta_12 <- 1*(t_12_n <= t_1D_n) #*delta_01
  delta_1D <- 1-delta_12
  delta_12[delta_01==0] <- 0
  delta_1D[delta_01==0] <- 0
  
  ### Death time
  t_D_n <- t_01_n + t_12_n + t_2D_n
  t_D_n[delta_0D==1] <- t_0D_n[delta_0D==1]
  t_D_n[delta_1D==1] <- (t_01_n + t_1D_n)[delta_1D==1]
  
  
  ### Ill time
  t_1_n <- t_01_n 
  t_1_n[delta_01==0] <- Inf
  
  ### Disease progression time
  t_2_n <- t_01_n + t_12_n
  t_2_n[delta_12==0] <- Inf
  
  res.list <- list(illtime1=t_1_n, illtime2=t_2_n, deathtime=t_D_n)
  return(res.list)
}

#illtime1 <- t_1_n
#illtime2 <- t_2_n
#lifetime <- t_D_n 
#b0 <- 0; b1 <- 10; S0 <- b1

# generate a simple random sample "population" 
gen_popu.f <- function(n, X, b0, b1, S0, illtime1, illtime2, lifetime){
  
  ### generate birth date
  B <- runif(n, min=b0, max=b1)
  A0 <- S0 - B # recruitment age
  
  ### observation of 0 to 1 transition
  delta1 <- 1*(illtime1<illtime2)*(illtime1<lifetime) #1*(illtime < lifetime)
  delta2 <- 1*(illtime1<illtime2)*(illtime2<lifetime)
  
  ### Y \in {0, 1}
  Y1 <- 1*(illtime1 < A0)*(lifetime > A0)
  Y2 <- 1*(illtime2 < A0)*(lifetime > A0)
  
  id <- 1:n
  
  res.df <- cbind(id, X, A0, Y1, Y2, illtime1, illtime2, lifetime, delta1, delta2)
  colnames(res.df) <- c("id", "X1", "X2", "A0", "Y1", "Y2", "A1", "A2", "AD", "delta1", "delta2")
  res.df <- as.data.frame(res.df)
  
  return(res.df)
}


# test
#adcen <- Inf
#rcen.param <- 0.05 #labd_CR
##res_popu <- res.df
#dt_Y1 = subset(res_popu, Y1==1 & Y2==0)
#dt_Y2 = subset(res_popu, Y2==1)
#n1_samp_Y1 = n1_samp_Y2 = 1000
#R_Y1 <- rep(0, n1_samp_Y1); R_Y2 = rep(0, n1_samp_Y2)
#R_Y1 [sample(1:n1_samp_Y1, 500)] <- 1
#R_Y2[sample(1:n1_samp_Y2, 500)] <- 1

# select a phase I sample
gen_df_obs.f <- function(n1_samp_Y1, n1_samp_Y2, R_Y1, R_Y2, adcen, rcen.param, dt_Y1, dt_Y2){
  
  ### generate random censoring time from exponential dist.
  if(is.finite(rcen.param)){
    random_CR_Y1 <- rexp(n1_samp_Y1, rate=rcen.param)
    random_CR_Y2 <- rexp(n1_samp_Y2, rate=rcen.param)
  }else{
    random_CR_Y1 = adcen +1e-04
    random_CR_Y2 = adcen +1e-04
  }
  
  s1_id_Y1 <- sample(1:length(dt_Y1$id), n1_samp_Y1)
  s1_id_Y2 <- sample(1:length(dt_Y2$id), n1_samp_Y2) #+ n1_samp_Y1
  
  ### observed censoring time
  A0_Y1 <- dt_Y1$A0[s1_id_Y1]
  A0_Y2 <- dt_Y2$A0[s1_id_Y2]
  obs_CR_Y1 <- pmin(adcen, random_CR_Y1) + A0_Y1
  obs_CR_Y2 <- pmin(adcen, random_CR_Y2) + A0_Y2
  
  AD_Y1 = dt_Y1$AD[s1_id_Y1]; AD_Y2 = dt_Y2$AD[s1_id_Y2]
  A1_Y1 = dt_Y1$A1[s1_id_Y1]; A1_Y2 = dt_Y2$A1[s1_id_Y2]
  A2_Y1 = dt_Y1$A2[s1_id_Y1]; A2_Y2 = dt_Y2$A2[s1_id_Y2]
  
  ### observed death time
  A_dgr_Y1 <- pmin(AD_Y1, obs_CR_Y1); A_dgr_Y2 <- pmin(AD_Y2, obs_CR_Y2)

  ### observed illtime1
  A1_dgr_Y1 <- pmin(A1_Y1, A_dgr_Y1); A1_dgr_Y2 <- pmin(A1_Y2, A_dgr_Y2)
  
  ### observed illtime2
  A2_dgr_Y1 <- pmin(A2_Y1, A_dgr_Y1); ############## not correct 
  A2_dgr_Y2 <- pmin(A2_Y2, A_dgr_Y2)
  
  
  T0_Y1 = A0_Y1 - A1_Y1; T0_Y2 = A0_Y2 - A1_Y2
  T2_dgr_Y1 = A2_dgr_Y1 - A1_Y1; T2_dgr_Y2 = A2_dgr_Y2 - A1_Y2 
  T_dgr_Y1 = A_dgr_Y1 - A1_Y1; T_dgr_Y2 = A_dgr_Y2 - A1_Y2 
  
  ### sojourn time in state 2
  #T2_dgr_Y1 <- A2_dgr_Y1 - A1_dgr_Y1
  
  ### observation of fatal transition
  delta_D_n_Y1 <- 1*(AD_Y1 <= obs_CR_Y1); delta_D_n_Y2 <- 1*(AD_Y2 <= obs_CR_Y2)
  
  ### observation of 0 to 1 transition
  delta_1_n_Y1 <- 1*(A1_Y1 <= A1_dgr_Y1); delta_1_n_Y2 <- 1*(A1_Y2 <= A1_dgr_Y2)
  
  ### observation of 1 to 2 transition
  delta_2_n_Y1 <- 1*(A2_Y1 <= A2_dgr_Y1); delta_2_n_Y2 <- 1*(A2_Y2 <= A2_dgr_Y2)
  
  id_Y1 <- 1:n1_samp_Y1; id_Y2 <- (1:n1_samp_Y2) + n1_samp_Y1
  
  X1_Y1 <- dt_Y1$X1[s1_id_Y1]; X1_Y2 <- dt_Y2$X1[s1_id_Y2]
  X1_Y1[R_Y1==0] <- NA; X1_Y2[R_Y2==0] <- NA
  X2_Y1 <- dt_Y1$X2[s1_id_Y1]; X2_Y2 <- dt_Y2$X2[s1_id_Y2]
  
  Y1 = ifelse(A0_Y1 > A1_Y1 & A0_Y1 < A2_Y1, 1, 2)
  Y2 = ifelse(A0_Y2 > A1_Y2 & A0_Y2 > A2_Y2, 2, 1)
  
  res.df_Y1 <- cbind(id_Y1, R_Y1, Y1, A0_Y1, T0_Y1, A1_dgr_Y1, A2_dgr_Y1, T2_dgr_Y1, A_dgr_Y1, T_dgr_Y1, delta_1_n_Y1, delta_2_n_Y1, delta_D_n_Y1, X1_Y1, X2_Y1)
  res.df_Y2 <- cbind(id_Y2, R_Y2, Y2, A0_Y2, T0_Y2, A1_dgr_Y2, A2_dgr_Y2, T2_dgr_Y2, A_dgr_Y2, T_dgr_Y2, delta_1_n_Y2, delta_2_n_Y2, delta_D_n_Y2, X1_Y2, X2_Y2)
  
  res.df <- rbind(res.df_Y1, res.df_Y2)
  
  colnames(res.df) <- c("id", "R", "Y", "A0", "T0", "A1_dagger", "A2_dagger", "T2_dagger", "A_dagger","T_dagger", "delta1", "delta2", "deltaD", "X1", "X2")
  res.df <- as.data.frame(res.df)
  
  #res.df.list = list(df_Y1=res.df_Y1,
  #                   df_Y2=res.df_Y2)
  
  return(res.df)
}


#head(do.call("cbind", res.list))
#sum(is.finite(res.list$illtime1))
#sum(is.finite(res.list$illtime2))
#mean(res.list$illtime1[is.finite(res.list$illtime1)])
#mean(res.list$illtime2[is.finite(res.list$illtime2)])
#mean(res.list$deathtime[is.finite(res.list$deathtime)])




