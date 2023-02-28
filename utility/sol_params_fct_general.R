
# the joint distrition for (X1, X2)
sol.probvec <- function(P1, P2, OR)
{
  # find p11
  g = function(p) log(p) + log(1-P1-P2+p) - log(P1-p) - log(P2-p) - log(OR)
  br = c( max(0,P1+P2-1), min(P1,P2) ) 
  p11 = uniroot(g, br)$root
  
  # fill in other cell probabilities
  p10 = P1 - p11
  p01 = P2 - p11
  p00 = 1-p11-p10-p01
  
  prob.vec <- c(p11,p10,p01,p00)
  
  return(prob.vec)
}

hazard_T.f <- function(t, X, coeffs.beta, family="weibull", param){
  
  if(family=="weibull"){
    
    scale_pm <- param[1]
    h.t <- scale_pm*param[2]*(scale_pm*t)^(param[2]-1)*exp(as.vector(X%*% coeffs.beta))
    #exp(sum(X*coeffs.beta))*param[1]*param[2]*(param[1]*t)^(param[2]-1)
  }
  
  if(family=="pwc"){
    
    dim_pm <- length(param)
    brks_vec <- param[-(1:((dim_pm+1)/2))]
    pwc_vec <- param[1:((dim_pm+1)/2)]
    
    h.t <- hpch(t, cuts=brks_vec,levels=pwc_vec)*exp(as.vector(X%*% coeffs.beta))
    #f.t <- dpch(t, cuts=brks_vec, levels=pwc_vec)
    #survF.t <- 1-ppch(t, cuts=brks_vec, levels=pwc_vec)
    #h.t <- f.t/survF.t
  }
  
  return(h.t)
  
}

surv_T.f <- function(t, X, coeffs.beta, family="weibull", param){
  
  if(family=="weibull"){
    
    scale_pm <- param[1]
    shape_pm <- param[2]
    
    survF <- exp(- (scale_pm*t)^(shape_pm)*exp(as.vector(X%*% coeffs.beta)))
    
  }
  
  if(family=="pwc"){
    dim_pm <- length(param)
    brks_vec <- param[-(1:((dim_pm+1)/2))]
    pwc_vec <- param[1:((dim_pm+1)/2)]
    
    survF <- exp(-Hpch(t, cuts=brks_vec, levels=pwc_vec)*exp(as.vector(X %*% coeffs.beta)))   #1- ppch(t, cuts=brks_vec, levels=pwc_vec)^(exp(as.vector(X %*% coeffs.beta)))
  }
  
  return(survF)
  
}

# gaussian legendre 
# prepare gaussleg.f(): reference of Jooyoung's code
gaussleg.f <- function(n, x1, x2) {
  EPS <- 3e-14
  
  m <- (n + 1)/2
  xm <- 0.5 * (x2 + x1)
  xl <- 0.5 * (x2 - x1)
  
  x <- rep(0, n)
  w <- rep(0, n)
  for (i in 1:m) {
    z <- cos(pi * (i - 0.25)/(n + 0.5))
    
    tol <- 9999
    while (tol > EPS) {
      p1 <- 1
      p2 <- 0
      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- ((2 * j - 1) * z * p2 - (j - 1) * p3)/j
      }
      
      pp <- n * (z * p1 - p2)/(z * z - 1)
      z1 <- z
      z <- z1 - p1/pp
      
      tol <- abs(z - z1)
      if (tol <= EPS) {
        break
      }
    }
    
    x[i] = xm - xl * z
    x[n + 1 - i] = xm + xl * z
    w[i] = (2 * xl)/((1 - z * z) * pp * pp)
    w[n + 1 - i] = w[i]
  }
  
  return(as.matrix(data.frame(location = x, weight = w)))
}


# constraints: P(N_01(inf)=1)=0.1 or 0.25
#              P(N_cdot D(inf)=1)=0.99

P_00_g_X.f <- function(t, X, coeffs.beta_T1, coeffs.beta_TD, 
                       family_T1="weibull", param_T1, 
                       family_TD="weibull", param_TD){
  
  survF_T1 = surv_T.f(t,X,coeffs.beta_T1, family_T1, param_T1)
  survF_T0D = surv_T.f(t,X,coeffs.beta_TD, family_TD, param_TD)
  
  res <- survF_T1*survF_T0D
  return(res)
}

P_11_g_X.f <- function(t1, t2, X, A1, coeffs.beta_T2, coeffs.beta_TD, gam, 
                       family_T2="weibull", param_T2,
                       family_TD="weibull", param_TD){
  
  survF_TD = surv_T.f(t2,X,coeffs.beta_TD, family_TD, param_TD)/surv_T.f(t1,X,coeffs.beta_TD, family_TD, param_TD)
  survF_T2 = surv_T.f(t2-t1, X, coeffs.beta_T2[-1], family_T2, param_T2)
  res <- survF_TD^(exp(gam))*survF_T2^(exp(coeffs.beta_T2[1]*log(A1)))
  return(res)
  
}

P_22_g_X.f <- function(t1, t2, X, coeffs.beta_TD, gam, family_TD="weibull", param_TD){
  
  survF_TD = surv_T.f(t2,X,coeffs.beta_TD, family_TD, param_TD)/surv_T.f(t1,X,coeffs.beta_TD, family_TD, param_TD)
  
  res <- survF_TD^(exp(gam))
  return(res)
  
}


P_Y0_g_X.f <- function(X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_TD, 
                              family_T1="weibull", param_T1, 
                              family_TD="weibull", param_TD){
  
  
  inner_g <- function(b){
    
    a = set_S0 - b
    P_00_g_X.f(a, X, coeffs.beta_T1, coeffs.beta_TD, 
               family_T1, param_T1, 
               family_TD, param_TD)
  }
  
  res <- integrate(inner_g, b0, b1)$val/(b1-b0)
  
  return(res)
}


P_Y1_g_X_gauleg.f <- function(X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                              family_T1="weibull", param_T1, 
                              family_T2="weibull", param_T2, 
                              family_TD="weibull", param_TD, gam){
  
  inner_g <- function(b, e){
    
    a0 = set_S0 - b
    a1 = e - b
    
    p00 <- P_00_g_X.f(a1, X, coeffs.beta_T1, coeffs.beta_TD, 
                      family_T1, param_T1, 
                      family_TD, param_TD)
    lam <- hazard_T.f(a1, X, coeffs.beta_T1, family_T1, param_T1)
    p11 <- P_11_g_X.f(a1, a0, X, a1, coeffs.beta_T2, coeffs.beta_TD, gam, family_T2, param_T2, family_TD, param_TD)
    p00*lam*p11
  }
  
  gau_b = gaussleg.f(20, b0, b1)
  u_b = gau_b[,1]; w_b = gau_b[,2]
  
  uu_e = NULL
  ww_e = NULL
  for(k in 1:20){
    
    b = u_b[k]
    
    gau_e = gaussleg.f(20, b, set_S0)
    u_e = gau_e[,1]; w_e = gau_e[,2]
    
    uu_e = c(uu_e, u_e)
    ww_e = c(ww_e, w_e)
  }
  
  uu_b = rep(u_b, each=20); ww_b=rep(w_b, each=20)
  
  res = sum(inner_g(uu_b, uu_e)*ww_b*ww_e)/(b1-b0)
  
  return(res)
}


P_Y2_g_X_gauleg.f <- function(X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2,coeffs.beta_TD, 
                              family_T1="weibull", param_T1, 
                              family_T2="weibull", param_T2, 
                              family_TD="weibull", param_TD, gam1, gam2){
  
  inner_g2 <- function(b, e1, e2){
  #inner_g2 <- function(a0, a1, a2){  
    a0 = set_S0 - b
    a1 = e1 - b
    a2 = e2 - b
    
    p00 <- P_00_g_X.f(a1, X, coeffs.beta_T1, coeffs.beta_TD, 
                      family_T1, param_T1, 
                      family_TD, param_TD)
    lam <- hazard_T.f(a1, X, coeffs.beta_T1, family_T1, param_T1)
    p11 <- P_11_g_X.f(a1, a2, X, a1, coeffs.beta_T2, coeffs.beta_TD, gam1, family_T2, param_T2, family_TD, param_TD)
    h <- hazard_T.f(a2-a1, X, coeffs.beta_T2[-1], family_T2, param_T2)*exp(coeffs.beta_T2[1]*log(a1))
    p22 <- P_22_g_X.f(a2,a0, X, coeffs.beta_TD, gam2, family_TD, param_TD)
    p00*lam*p11*h*p22
  }
  
  inner_g <- function(b){
    
    gau_e1 = gaussleg.f(20, b, set_S0)
    u_e1 = gau_e1[,1]; w_e1 = gau_e1[,2]
    
    uu_e2 = NULL
    ww_e2 = NULL
    for(k in 1:20){
      
      e1 = u_e1[k]
      
      gau_e2 = gaussleg.f(20, e1, set_S0)
      u_e2 = gau_e2[,1]; w_e2 = gau_e2[,2]
      
      uu_e2 = c(uu_e2, u_e2)
      ww_e2 = c(ww_e2, w_e2)
    }
    
    uu_e1 = rep(u_e1, each=20); ww_e1=rep(w_e1, each=20)
    
    sum(inner_g2(b,uu_e1, uu_e2)*ww_e1*ww_e2)
  }
  
  gau_b = gaussleg.f(20, b0, b1)
  u_b = gau_b[,1]; w_b = gau_b[,2]
  #res <- sum(inner_g(u_b)*w_b)/(b1-b0)
  res <- sum(sapply(u_b, inner_g)*w_b)/(b1-b0)
  
  return(res)
}


solve_lam0_alp0_mu0_gauleg.f <- function(probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                         family_T1, family_T2, family_TD, 
                                         kai1, kai2, kaiD, gamma1, gamma2, 
                                         pY0, pY1, pY2, set_S0){
  
  obj.f <- function(y){
    
    x <- exp(y)
    
    P_Y0 = (probvec[1]*P_Y0_g_X.f(c(1,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_TD, 
                                  family_T1, param_T1=c(x[1], kai1), 
                                  family_TD, param_TD=c(x[3], kaiD))
            + probvec[2]*P_Y0_g_X.f(c(1,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_TD,
                                    family_T1, param_T1=c(x[1], kai1), 
                                    family_TD, param_TD=c(x[3], kaiD))
            + probvec[3]*P_Y0_g_X.f(c(0,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_TD,
                                    family_T1, param_T1=c(x[1], kai1), 
                                    family_TD, param_TD=c(x[3], kaiD))
            + probvec[4]*P_Y0_g_X.f(c(0,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_TD,
                                    family_T1, param_T1=c(x[1], kai1), 
                                    family_TD, param_TD=c(x[3], kaiD)))
    
    P_Y1 = (probvec[1]*P_Y1_g_X_gauleg.f(c(1,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                         family_T1, param_T1=c(x[1], kai1), 
                                         family_T2, param_T2=c(x[2], kai2), 
                                         family_TD, param_TD=c(x[3], kaiD),
                                         gamma1)
            + probvec[2]*P_Y1_g_X_gauleg.f(c(1,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                           family_T1, param_T1=c(x[1], kai1), 
                                           family_T2, param_T2=c(x[2], kai2), 
                                           family_TD, param_TD=c(x[3], kaiD),
                                           gamma1)
            + probvec[3]*P_Y1_g_X_gauleg.f(c(0,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                           family_T1, param_T1=c(x[1], kai1), 
                                           family_T2, param_T2=c(x[2], kai2), 
                                           family_TD, param_TD=c(x[3], kaiD),
                                           gamma1)
            + probvec[4]*P_Y1_g_X_gauleg.f(c(0,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                           family_T1, param_T1=c(x[1], kai1), 
                                           family_T2, param_T2=c(x[2], kai2), 
                                           family_TD, param_TD=c(x[3], kaiD),
                                           gamma1))
    
    P_Y2 = (probvec[1]*P_Y2_g_X_gauleg.f(c(1,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                      family_T1, param_T1=c(x[1], kai1), 
                                      family_T2, param_T2=c(x[2], kai2), 
                                      family_TD, param_TD=c(x[3], kaiD),
                                      gamma1, gamma2)
            + probvec[2]*P_Y2_g_X_gauleg.f(c(1,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                        family_T1, param_T1=c(x[1], kai1), 
                                        family_T2, param_T2=c(x[2], kai2), 
                                        family_TD, param_TD=c(x[3], kaiD),
                                        gamma1, gamma2)
            + probvec[3]*P_Y2_g_X_gauleg.f(c(0,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                        family_T1, param_T1=c(x[1], kai1), 
                                        family_T2, param_T2=c(x[2], kai2), 
                                        family_TD, param_TD=c(x[3], kaiD),
                                        gamma1, gamma2)
            + probvec[4]*P_Y2_g_X_gauleg.f(c(0,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                        family_T1, param_T1=c(x[1], kai1), 
                                        family_T2, param_T2=c(x[2], kai2), 
                                        family_TD, param_TD=c(x[3], kaiD),
                                        gamma1, gamma2) )
    
    #print(c(res1, res2))
    res0 = pY0 - P_Y0
    res1 = pY1 - P_Y1
    res2 = pY2 - P_Y2
    
    #res <- crossprod(c(res0, res1, res2))#sqrt(abs(res1)+abs(res2))
    res <- c(res0, res1, res2)
    #print(res)
    return(res)
  }
  
  
  #res <- optim(c(log(1),log(1)), obj.f, method="L-BFGS-B")$par
  res <- nleqslv(rep(0.01, 3), obj.f, 
                 control = list(ftol=1e-09, xtol=1e-09), method="Newton")$x
  
  
  return(res)
  
}

P_Y1_0D_g_X_gauleg.f <- function(X, C, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                              family_T1="weibull", param_T1, 
                              family_T2="weibull", param_T2, 
                              family_TD="weibull", param_TD, gam){
  
  inner_g2 <- function(b, e, cd){
    
    a0 = set_S0 - b
    a1 = e - b
    acd = a0 + cd
    
    p00 <- P_00_g_X.f(a1, X, coeffs.beta_T1, coeffs.beta_TD, 
                      family_T1, param_T1, 
                      family_TD, param_TD)
    lam <- hazard_T.f(a1, X, coeffs.beta_T1, family_T1, param_T1)
    p11 <- P_11_g_X.f(a1, acd, X, a1, coeffs.beta_T2, coeffs.beta_TD, gam, family_T2, param_T2, family_TD, param_TD)
    eta <- hazard_T.f(acd, X, coeffs.beta_TD, family_TD, param_TD)*exp(gam)
    
    p00*lam*p11*eta
  }
  
  inner_g <- function(cd){
    
    gau_b = gaussleg.f(20, b0, b1)
    u_b = gau_b[,1]; w_b = gau_b[,2]
    
    uu_e = NULL
    ww_e = NULL
    for(k in 1:20){
      
      b = u_b[k]
      
      gau_e = gaussleg.f(20, b, set_S0)
      u_e = gau_e[,1]; w_e = gau_e[,2]
      
      uu_e = c(uu_e, u_e)
      ww_e = c(ww_e, w_e)
    }
    
    uu_b = rep(u_b, each=20); ww_b=rep(w_b, each=20)
    
    res = sum(inner_g2(uu_b, uu_e, cd)*ww_b*ww_e)#/(b1-b0)
    
    res
  }
  
  gau_c = gaussleg.f(20, 0, C)
  u_c = gau_c[,1]; w_c = gau_c[,2]
  #res <- sum(inner_g(u_b)*w_b)/(b1-b0)
  res <- sum(sapply(u_c, inner_g)*w_c)/(b1-b0)
  
  return(res)
}

P_Y1_01_g_X_gauleg.f <- function(X, C, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                 family_T1="weibull", param_T1, 
                                 family_T2="weibull", param_T2, 
                                 family_TD="weibull", param_TD, gam){
  
  inner_g2 <- function(b, e){
    
    a0 = set_S0 - b
    a1 = e - b
    ac = a0 + C
    
    p00 <- P_00_g_X.f(a1, X, coeffs.beta_T1, coeffs.beta_TD, 
                      family_T1, param_T1, 
                      family_TD, param_TD)
    lam <- hazard_T.f(a1, X, coeffs.beta_T1, family_T1, param_T1)
    p11 <- P_11_g_X.f(a1, ac, X, a1, coeffs.beta_T2, coeffs.beta_TD, gam, family_T2, param_T2, family_TD, param_TD)
    #eta <- hazard_T.f(acd, X, coeffs.beta_TD, family_TD, param_TD)*exp(gam)
    
    p00*lam*p11
  }
  
  gau_b = gaussleg.f(20, b0, b1)
  u_b = gau_b[,1]; w_b = gau_b[,2]
  
  uu_e = NULL
  ww_e = NULL
  for(k in 1:20){
    
    b = u_b[k]
    
    gau_e = gaussleg.f(20, b, set_S0)
    u_e = gau_e[,1]; w_e = gau_e[,2]
    
    uu_e = c(uu_e, u_e)
    ww_e = c(ww_e, w_e)
  }
  
  uu_b = rep(u_b, each=20); ww_b=rep(w_b, each=20)
  
  res = sum(inner_g2(uu_b, uu_e)*ww_b*ww_e)/(b1-b0)
  
  return(res)
}


P_del2_g_Y1_X_gauleg.f <- function(set_CA, X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                   family_T1, param_T1, 
                                   family_T2, param_T2, 
                                   family_TD, param_TD, gam1){
  
  res1 = P_Y1_01_g_X_gauleg.f(X, set_CA, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                              family_T1, param_T1, 
                              family_T2, param_T2, 
                              family_TD, param_TD, gam1)
  
  res2 = P_Y1_0D_g_X_gauleg.f(X, set_CA, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                              family_T1, param_T1, 
                              family_T2, param_T2, 
                              family_TD, param_TD, gam1)
  
  res0 = P_Y1_g_X_gauleg.f(X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                           family_T1, param_T1, 
                           family_T2, param_T2, 
                           family_TD, param_TD, gam1)
  res = (res1 + res2)/res0
  
  return(1-res)
}

P_del2_g_Y1_X_gauleg_rcen.f <- function(lamC, X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                   family_T1, param_T1, 
                                   family_T2, param_T2, 
                                   family_TD, param_TD, gam1){
  
  gau_c = gaussleg.f(20, 0, 1)
  u_c = gau_c[,1]; w_c = gau_c[,2]
  
  inner_g1 <- function(c){
    
    res1 = P_Y1_01_g_X_gauleg.f(X, c, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                family_T1, param_T1, 
                                family_T2, param_T2, 
                                family_TD, param_TD, gam1)
    
    res2 = P_Y1_0D_g_X_gauleg.f(X, c, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                family_T1, param_T1, 
                                family_T2, param_T2, 
                                family_TD, param_TD, gam1)
    
    res = (res1 + res2)*lamC*exp(-lamC*c)
    
  }
  
  inner_g2 <- function(c){
    
    res1 = P_Y1_01_g_X_gauleg.f(X, 1/c, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                family_T1, param_T1, 
                                family_T2, param_T2, 
                                family_TD, param_TD, gam1)
    
    res2 = P_Y1_0D_g_X_gauleg.f(X, 1/c, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                family_T1, param_T1, 
                                family_T2, param_T2, 
                                family_TD, param_TD, gam1)
    
    res = (res1 + res2)*lamC*exp(-lamC*1/c)*(1/c)*(1/c)
    
  }
  
  res1 <- sum(sapply(u_c, inner_g1)*w_c)
  res2 <- sum(sapply(u_c, inner_g2)*w_c)
  
  res0 = P_Y1_g_X_gauleg.f(X, set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                           family_T1, param_T1, 
                           family_T2, param_T2, 
                           family_TD, param_TD, gam1)
  
  res = (res1+res2)/res0
  
  return(1-res)
}

solve_set_CA_gauleg.f <- function(q_C, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                  family_T1="weibull", param_T1, 
                                  family_T2="weibull", param_T2, 
                                  family_TD="weibull", param_TD, gam1){
  
  
  obj.f <- function(x){
    
    CA <- exp(x)
    
    P02_11 <- P_del2_g_Y1_X_gauleg.f(CA, c(1,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    P02_10 <- P_del2_g_Y1_X_gauleg.f(CA, c(1,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    P02_01 <- P_del2_g_Y1_X_gauleg.f(CA, c(0,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    P02_00 <- P_del2_g_Y1_X_gauleg.f(CA, c(0,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    res <- (probvec[1]*(P02_11)
            +probvec[2]*(P02_10)
            +probvec[3]*(P02_01)
            +probvec[4]*(P02_00)) - q_C
    
    #print(res)
    #print(CA)
    
    return((res))
  }
  
  
  res <- uniroot(obj.f, c(-100, 100))$root
  #res <- nleqslv(log(0.01), obj.f, control = list(ftol=1e-09, xtol=1e-09), method="Newton")$x
  #res <- optim(log(0.01), obj.f, method="L-BFGS-B")$par
  
  return(res)
  
}

solve_lamC_gauleg.f <- function(q_C, set_S0, probvec, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD,
                                  family_T1="weibull", param_T1, 
                                  family_T2="weibull", param_T2, 
                                  family_TD="weibull", param_TD, gam1){
  
  
  obj.f <- function(x){
    
    lamC <- exp(x)
    
    P02_11 <- P_del2_g_Y1_X_gauleg_rcen.f(lamC, c(1,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    P02_10 <- P_del2_g_Y1_X_gauleg_rcen.f(lamC, c(1,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    P02_01 <- P_del2_g_Y1_X_gauleg_rcen.f(lamC, c(0,1), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    P02_00 <- P_del2_g_Y1_X_gauleg_rcen.f(lamC, c(0,0), set_S0, b0, b1, coeffs.beta_T1, coeffs.beta_T2, coeffs.beta_TD, 
                                     family_T1, param_T1, 
                                     family_T2, param_T2, 
                                     family_TD, param_TD, gam1)
    res <- (probvec[1]*(P02_11)
            +probvec[2]*(P02_10)
            +probvec[3]*(P02_01)
            +probvec[4]*(P02_00)) - q_C
    
    #print(res)
    #print(CA)
    
    return((res))
  }
  
  
  res <- uniroot(obj.f, c(-10, 10))$root
  #res <- nleqslv(log(0.01), obj.f, control = list(ftol=1e-09, xtol=1e-09), method="Newton")$x
  #res <- optim(log(0.01), obj.f, method="L-BFGS-B")$par
  
  return(res)
  
}
