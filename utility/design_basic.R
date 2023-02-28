
expit <- function(x){
  exp(x)/(1+exp(x))
}

phase_I_stratify_IDM.f <- function(pcuts=NULL, dt, design.factor){
  
  style = paste0(design.factor, collapse="_")
  
  if(is.null(pcuts)){
    pcuts <- c(0.33, 0.66)
  }
  
  T0_cut <- quantile(dt$T0, p=pcuts)
  T0bar <- cut(dt$T0, c(-0.01,T0_cut,Inf), labels = paste0(1:(length(T0_cut)+1)))
  
  A1_cut <- quantile(dt$A1_dagger, p=pcuts)
  A1bar <- cut(dt$A1_dagger, c(-0.01,A1_cut,Inf), labels = paste0(1:(length(A1_cut)+1)))
  
  T2dgr_cut <- quantile(dt$T2_dagger, p=pcuts)
  T2dgrbar <- cut(dt$T2_dagger, c(-0.01,T2dgr_cut,Inf), labels = paste0(1:(length(T2dgr_cut)+1)))
  
  dt$T0bar <- T0bar
  dt$A1bar <- A1bar
  dt$T2dgrbar <- T2dgrbar
  
  if(style=="T0"){
    z <- xtabs(~T0bar, dt)
    group <- factor(paste0(dt$T0bar))
    ref.group <- 1:(length(pcuts)+1)
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="A1"){
    
    z <- xtabs(~A1bar, dt)
    group <- factor(paste0(dt$A1bar))
    ref.group <- 1:(length(pcuts)+1)
    strata.n <- as.data.frame(z)$Freq # give size of each strata
  }
  
  if(style=="A1_X2"){
    z <- xtabs(~A1bar+X2, dt)
    group <- factor(paste0(dt$A1bar, dt$X2))
    group.mat <- expand.grid(1:3, 0:1)
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
  }
  
  if(style=="del2_A1_X2"){
    z <- xtabs(~delta2+A1bar+X2, dt)
    group <- factor(paste0(dt$delta2,dt$A1bar, dt$X2))
    group.mat <- expand.grid(0:1, 1:3, 0:1)
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
  }
  
  if(style=="Y"){
    z <- xtabs(~Y, dt)
    group <- factor(paste0(dt$Y))
    ref.group <- 1:2
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="del2"){
    z <- xtabs(~delta2, dt)
    group <- factor(paste0(dt$delta2))
    ref.group <- 0:1
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="Y_del2"){
    z <- xtabs(~Y+delta2, dt)
    group <- factor(paste0(dt$Y,dt$delta2))
    group.mat <- expand.grid(1:2, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="Y_del2_X2"){
    z <- xtabs(~Y+delta2+X2, dt)
    group <- factor(paste0(dt$Y,dt$delta2,dt$X2))
    group.mat <- expand.grid(1:2, 0:1, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="Y_del2_delD"){
    z <- xtabs(~Y+delta2+deltaD, dt)
    group <- factor(paste0(dt$Y,dt$delta2, dt$deltaD))
    group.mat <- expand.grid(1:2, 0:1, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], 
                        group.mat[,3])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="del2_T2dgr"){
    z <- xtabs(~delta2+T2dgrbar, dt)
    group <- factor(paste0(dt$delta2, dt$T2dgrbar))
    group.mat <- expand.grid(0:1, 1:3)
    #ind = (group.mat[,1]==2 & group.mat[,2]==0)
    #group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
  }
  
  if(style=="del2_T0"){
    z <- xtabs(~delta2+T0bar, dt)
    group <- factor(paste0(dt$delta2, dt$T0bar))
    group.mat <- expand.grid(0:1, 1:3)
    #ind = (group.mat[,1]==2 & group.mat[,2]==0)
    #group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
  }

  if(style=="Y_T0"){
    z <- xtabs(~Y+T0bar, dt)
    group <- factor(paste0(dt$Y, dt$T0bar))
    group.mat <- expand.grid(1:2, 1:(length(pcuts)+1))
    #ind = (group.mat[,1]==2 & group.mat[,2]==0)
    #group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="Y_T2dgr"){
    z <- xtabs(~Y+T2dgrbar, dt)
    group <- factor(paste0(dt$Y, dt$T2dgrbar))
    group.mat <- expand.grid(1:2, 1:(length(pcuts)+1))
    #ind = (group.mat[,1]==2 & group.mat[,2]==0)
    #group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="Y_del2_T0"){
    z <- xtabs(~Y+delta2+T0bar, dt)
    group <- factor(paste0(dt$Y, dt$delta2, dt$T0bar))
    group.mat <- expand.grid(1:2, 0:1, 1:3)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="Y_del2_A1"){
    z <- xtabs(~Y+delta2+A1bar, dt)
    group <- factor(paste0(dt$Y, dt$delta2, dt$A1bar))
    group.mat <- expand.grid(1:2, 0:1, 1:3)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="Y_del2_T2dgr"){
    z <- xtabs(~Y+delta2+T2dgrbar, dt)
    group <- factor(paste0(dt$Y, dt$delta2, dt$T2dgrbar))
    group.mat <- expand.grid(1:2, 0:1, 1:3)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  } 
  
  if(style=="del2_delD_T2dgr"){
    z <- xtabs(~delta2+deltaD+T2dgrbar, dt)
    group <- factor(paste0(dt$delta2,dt$deltaD,dt$T2dgrbar))
    group.mat <- expand.grid(0:1, 0:1, 1:3)
    #ind = (group.mat[,1]==2 & group.mat[,2]==0)
    #group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="Y_del2_delD_T0"){
    z <- xtabs(~Y+delta2+deltaD+T0bar, dt)
    group <- factor(paste0(dt$Y,dt$delta2, dt$deltaD, dt$T0bar))
    group.mat <- expand.grid(1:2, 0:1, 0:1, 1:3)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], 
                        group.mat[,3], group.mat[,4])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="Y_del2_delD_T2dgr"){
    z <- xtabs(~Y+delta2+deltaD+T2dgrbar, dt)
    group <- factor(paste0(dt$Y,dt$delta2, dt$deltaD, dt$T2dgrbar))
    group.mat <- expand.grid(1:2, 0:1, 0:1, 1:3)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], 
                        group.mat[,3], group.mat[,4])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
  }
  
  if(style=="Y_del2_delD_T0_A1_T2dgr_X2"){
    z <- xtabs(~Y+delta2+deltaD+T0bar+A1bar+T2dgrbar+X2, dt)
    group <- factor(paste0(dt$Y,dt$delta2, dt$deltaD, dt$T0bar,dt$A1bar,dt$T2dgrbar, dt$X2))
    group.mat <- expand.grid(1:2, 0:1, 0:1, 1:3, 1:3, 1:3, 0:1)
    ind = (group.mat[,1]==2 & group.mat[,2]==0)
    group.mat <- group.mat[!ind,]
    ref.group <- paste0(group.mat[,1], group.mat[,2], group.mat[,3],
                        group.mat[,4], group.mat[,5], group.mat[,6], group.mat[,7])
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
  }
  
  dt$group <- group
  
  res <- list(strata.n=strata.n,
              dt_ext = dt,
              ref.strata = ref.group)
  
  return(res)
  
}

phase_I_stratify_Mjmu_IDM.f <- function(pcuts=NULL, dt, Mjmu){
  
  if(is.null(pcuts)){
    pcuts <- c(0.33, 0.66)
  }
  
  Mjmu_cut <- quantile(Mjmu, p=pcuts)
  Mjmubar <- cut(Mjmu, c(-Inf,Mjmu_cut,Inf), labels = paste0(1:(length(Mjmu_cut)+1)))
  
  dt$Mjmubar <- Mjmubar
  #z <- xtabs(~Y+Mjmubar, dt)
  #group <- factor(paste0(dt$Y,dt$Mjmubar))
  #group.mat <- expand.grid(1:2, 1:(length(pcuts)+1))
  
  #ref.group <- paste0(group.mat[,1], group.mat[,2])
  #group <- factor(paste0(dt$Y, dt$Mjmubar))
  
  z <- xtabs(~Mjmubar, dt)
  group <- factor(dt$Mjmubar)
  ref.group <- 1:(length(pcuts)+1)
  strata.n <- as.data.frame(z)$Freq # give size of each strata
  
  dt$group <- group
  
  res <- list(strata.n=strata.n,
              dt_ext = dt,
              ref.strata = ref.group)
  
  return(res)
  
}

subopt.f <- function(n2_samp, design="srs", phaseI_res){
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- length(dt_ext$id)
  
  if(design=="srs"){
    
    s_id <- as.numeric(sample(as.character(dt_ext$id),n2_samp))
    
  }
  
  if(design=="bal"){
    
    num_strata <- length(strata.n)
    n_circ_vec <- pmin(n2_samp/num_strata, strata.n)
    sum_n_circ <- sum(n_circ_vec)
    
    left <- n2_samp - sum_n_circ
    
    while(left >0 ){
      left_bal <- strata.n - n_circ_vec
      add_bal <- numeric(num_strata)
      for(i in 1:num_strata){
        div <- sum(left_bal >0)
        add_bal[i] <- ifelse(left_bal[i] < left/div,left_bal[i],left/div)
      }
      n_circ_vec <- n_circ_vec + add_bal
      left <- n2_samp - sum(n_circ_vec)
    }
    
    s_bal <- floor(n_circ_vec)
    
    s.prob <- s_bal/strata.n
    
    s_id <- unlist(sapply(which(s_bal!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group%in% ref.group[x]]),s_bal[x]))}))
    #if(sum(s_bal)<n2_samp){
    #  s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n2_samp-sum(s_bal))))
    #}
    
  }
  
  if(design=="ext_Y_6strata"){
    num_strata <- length(strata.n)
    # S_etr: 21, 13
    s_etr <- rep(0.05, num_strata)*strata.n
    etr.strata <- which(ref.group %in% c("21","13"))
    s_etr[etr.strata] <- floor(pmax(pmin(strata.n[etr.strata], (n2_samp-sum(s_etr))/length(etr.strata)+s_etr[etr.strata]),0))
    #s_etr[-etr.strata] <- round(pmin(strata.n[-etr.strata], propbal/6)) 
    naive_etr <- s_etr
    left <- n2_samp - sum(naive_etr)
    
    while(left >0 ){
      left_etr <- strata.n - naive_etr
      add_etr <- numeric(num_strata)
      for(i in 1:num_strata){
        div <- sum(left_etr >0)
        add_etr[i] <- ifelse(left_etr[i] < left/div,left_etr[i],left/div)
      }
      naive_etr <- naive_etr + add_etr
      left <- n2_samp - sum(naive_etr)
    }
    
    s_etr <- floor(naive_etr)
    
    #s.prob <- s_etr/strata.n
    
    s_id <- unlist(sapply(which(s_etr!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group%in% ref.group[x]]),s_etr[x]))}))
    if(sum(s_etr)<n2_samp){
      s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n2_samp-sum(s_etr))))
    }
  }
  
  if(design=="ext_del2_6strata"){
    num_strata <- length(strata.n)
    # S_etr: 21, 13
    s_etr <- rep(0.05, num_strata)*strata.n
    etr.strata <- which(ref.group %in% c("11","03"))
    s_etr[etr.strata] <- floor(pmax(pmin(strata.n[etr.strata], (n2_samp-sum(s_etr))/length(etr.strata)+s_etr[etr.strata]),0))
    #s_etr[-etr.strata] <- round(pmin(strata.n[-etr.strata], propbal/6)) 
    naive_etr <- s_etr
    left <- n2_samp - sum(naive_etr)
    
    while(left >0 ){
      left_etr <- strata.n - naive_etr
      add_etr <- numeric(num_strata)
      for(i in 1:num_strata){
        div <- sum(left_etr >0)
        add_etr[i] <- ifelse(left_etr[i] < left/div,left_etr[i],left/div)
      }
      naive_etr <- naive_etr + add_etr
      left <- n2_samp - sum(naive_etr)
    }
    
    s_etr <- floor(naive_etr)
    
    #s.prob <- s_etr/strata.n
    
    s_id <- unlist(sapply(which(s_etr!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group%in% ref.group[x]]),s_etr[x]))}))
    if(sum(s_etr)<n2_samp){
      s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n2_samp-sum(s_etr))))
    }
  }
  
  if(design=="ext_Y_T0"){
    
    T0 <- dt_ext$T0
    indY_1 <- which(dt_ext$Y==1); indY_2 <- which(dt_ext$Y==2)
    order_resi2 = order(T0[indY_2])
    order_resi1 = order(T0[indY_1])
    
    nY_2 <- length(indY_2); nY_1 <- length(indY_1)
    
    k = min(n2_samp/2, nY_1)
    #s_id_2 = dt$id[c(indY_1[order_resi1[(nY_1-k+1):nY_1]], indY_0[order_resi0[1:(n_sample-k)]])]
    s_id = dt_ext$id[c(indY_2[order_resi2[1:(n2_samp-k)]], indY_1[order_resi1[(nY_1 -k +1):nY_1]])]
    
  }
  
  if(design=="ext_del2_T0"){
    
    T0 <- dt_ext$T0
    indY_1 <- which(dt_ext$delta2==1); indY_2 <- which(dt_ext$delta2==0)
    order_resi2 = order(T0[indY_2])
    order_resi1 = order(T0[indY_1])
    
    nY_2 <- length(indY_2); nY_1 <- length(indY_1)
    
    k = min(n2_samp/2, nY_1)
    #s_id_2 = dt$id[c(indY_1[order_resi1[(nY_1-k+1):nY_1]], indY_0[order_resi0[1:(n_sample-k)]])]
    s_id = dt_ext$id[c(indY_2[order_resi2[1:(n2_samp-k)]], indY_1[order_resi1[(nY_1 -k +1):nY_1]])]
    
  }
  
  if(design=="ext_del2_T2dgr"){
    
    T0 <- dt_ext$T2dgr
    indY_1 <- which(dt_ext$delta2==1); indY_2 <- which(dt_ext$delta2==0)
    order_resi2 = order(T0[indY_2])
    order_resi1 = order(T0[indY_1])
    
    nY_2 <- length(indY_2); nY_1 <- length(indY_1)
    
    k = min(n2_samp/2, nY_1)
    #s_id_2 = dt$id[c(indY_1[order_resi1[(nY_1-k+1):nY_1]], indY_0[order_resi0[1:(n_sample-k)]])]
    s_id = dt_ext$id[c(indY_1[order_resi1[1:(n2_samp-k)]], indY_2[order_resi2[(nY_2 -k +1):nY_2]])]
    
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  #s.prob <- s_num/strata.n #
  
  if(design=="srs"){
    s.prob <- rep(n2_samp/sum(strata.n), length(strata.n))
  }else{
    s.prob <- s_num/strata.n #
  }  

  res <- list(s_id=s_id,
              s_prob=s.prob,
              s_m=s_num)
  
  return(res)
  
}

subopt_inf.f <- function(n2_samp, m, phaseI_res){
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- length(dt_ext$id)
  
  N2 = sum(dt_ext$group==2)
  n2 = m #*(1-s2_prop)
  
  N1 = sum(dt_ext$group==1)
  n1 = min(N1, (n2_samp-m)/2)
  n3 = n2_samp - n1 - n2
  
  s_id <- as.numeric(c(sample(as.character(subset(dt_ext, group==1)$id),n1),
                       sample(as.character(subset(dt_ext, group==2)$id),n2),
                       sample(as.character(subset(dt_ext, group==3)$id),n3)))
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n #
  
  res <- list(s_id=s_id,
              s_prob=s.prob,
              s_m=s_num)
  
  return(res)
}

scoreL_mu_noX.mortality.f <- function(dt){
  
  obj.f <- function(vartheta){
    
    params0_T2 = c(exp(vartheta[1]),1)
    params.reg_T2 = c(vartheta[1+1], 0, 0) #coeffs.beta_T1 #
    gam <- vartheta[2+1]
    eta = exp(vartheta[2+2])
    #params.X = vartheta[5+(1:2)] #c(eta0, eta1) #
    
    res1 = obslogLn_mortality_noX_vec(dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                                  "weibull", params0_T2, params.reg_T2, gam, eta)
    
    #res1 <- ifelse(is.finite(res1), res1, 0)
    
    res = -sum(res1)
    
    #print(res)
    #print(vartheta)
    #res0 <- scorePL.tr.f(gradstep=1e-06, vartheta, dt)
    
    #res <- colMeans(res0)
    return(res)
  }
  
  dim <- 4
  #est.res <- nlm(obj.f, p=rep(-0.1, dim), gradtol=1e-06, steptol=1e-06, hessian=TRUE)
  est.res <- optim(par=rep(-0.01, dim), fn=obj.f, method = c("L-BFGS-B"))
  #est.res <- nleqslv(rep(0.05, dim), obj.f, control = list(ftol=1e-09, xtol=1e-09), method="Newton", jacobian=TRUE)
  
  alp0 = exp(est.res$par[1])
  params.reg_T2 = c(est.res$par[1+1], 0, 0) #coeffs.beta_T1 #
  gam <- est.res$par[2+1]
  eta0 = exp(est.res$par[2+2])
  
  mu = params.reg_T2[1]*log(dt$A1_dagger) #+ params.reg[3]*dt_c$X2
  cov_reg = exp(mu)
  survF0 = exp(-alp0*dt$T0*cov_reg) #survF1_g_A1_X(T0, A1, 1, X2, family, params0, c(params_reg[1], 0, params_reg[3]))
  survG0 = exp(-eta0*dt$T0)
  h_star = (alp0*exp(mu)+eta0 - eta0*exp(gam))
  
  S1 = dt$delta2[dt$Y==1] - alp0*exp(mu[dt$Y==1])*( (dt$T2_dagger[dt$Y==1])  )
  S2 = 1 - alp0*exp(mu[dt$Y==2])*( (dt$T2_dagger[dt$Y==2])  )
  
  #I1 = (alp0*exp(mu)*dt$T2_dagger)[dt$Y==1]
  #I2 = (alp0*exp(mu)*dt$T2_dagger+ ((eta0*gam-eta0)/h_star)*alp0*exp(mu)/h_star)[dt$Y==2]
  
  aug_S1 = (alp0*exp(mu)*dt$T0)[dt$Y==1]
  aug_S2 =  (- exp(-h_star*dt$T0)/(1- exp(-h_star*dt$T0))*(alp0*exp(mu)*dt$T0))[dt$Y==2] + (alp0*exp(mu)/h_star)[dt$Y==2]-1
  
  #aug_I1 = (alp0*exp(mu)*dt$T0)[dt$Y==1]
  #ratio = exp(-h_star*dt$T0)/(1- exp(-h_star*dt$T0))
  #aug_I2 = (ratio*alp0*exp(mu)*dt$T0*(alp0*exp(mu)-ratio*alp0*exp(mu)*dt$T0))[dt$Y==2]
  
  #res = S1 + aug_0
  res = list(M1mu = S1, 
             M2mu = S2, 
             aug_S1 = aug_S1,
             aug_S2 = aug_S2
             #I1 = I1,
             #I2 = I2,
             #aug_I1 = aug_I1,
             #aug_I2 = aug_I2
             )
  return(res)
}

scorePL_mu_pwc.f <- function(dt, Xmat, main.cov, brks){
  
  if(any(is.finite(brks))){
    num.pwc = length(brks)+1
  }else{
    num.pwc=1
  }
  
  obs.f <- function(vartheta){
    
    params0_T2 = exp(vartheta[1:num.pwc])
    beta1 = 0#vartheta[num.pwc+1]
    params.reg_T2 = c(vartheta[num.pwc+(1:length(main.cov))])
    
    mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F] %*% params.reg_T2)
    
    res1 = obslogPLn_pwc_vec(dt$R, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$X1, brks, 
                             params0_T2, beta1, mu_X2_T2, rep(0, dim(dt)[1]))
    res = -sum(res1)#-sum(res1[is.finite(res1)]) #
    
    res <- ifelse(is.finite(res), res, NaN)
    #print(res)
    #print(vartheta)
    
    return(res)
  }
  
  dim <- num.pwc+length(main.cov)
  #est.res <- optim(par=rep(-0.01, dim), fn=obs.f, method = c("L-BFGS-B"), hessian = TRUE)
  est.res <- try(nlm(f=obs.f, p=rep(0.01, dim)))
  if(isTRUE(class(est.res)=="try-error")) {
    est.res <- try(nlm(f=obs.f, p=rep(-0.01, dim)))
    if(isTRUE(class(est.res)=="try-error")) {
      est.res <- try(nlm(f=obs.f, p=rep(0.01, dim)))
    }
  }
  vartheta = est.res$est

  params0_T2 = exp(vartheta[1:num.pwc])
  beta1 = 0#vartheta[num.pwc+1]
  params.reg_T2 = c(vartheta[num.pwc+(1:length(main.cov))])
  
  mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F] %*% params.reg_T2)
  
  logsurvF = survF1_g_A1_X2_pwc_vec(dt$T2_dagger, dt$A1_dagger, brks, params0_T2, mu_X2_T2,1)
  res = dt$delta2 + logsurvF
 
  return(res)
}

scoreL_mu_pwc_dm.f <- function(dt, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
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
    numpwc_T2D=1
  } 
  
  eta1=eta2=0; 
  mu_X2_T1D = rep(0, dim(Xmat)[1])
  mu_X2_T2D = rep(0, dim(Xmat)[1])
  mu_X = rep(0, dim(Xmat)[1])
  
  obs.f <- function(vartheta){
    
    params0_T2 = exp(vartheta[1:numpwc_T2])
    beta1 = 0#vartheta[numpwc_T2+1]
    params.reg_T2 = c(vartheta[numpwc_T2+(1:length(main.cov))])
    params0_T1D = exp(vartheta[numpwc_T2+length(main.cov)+(1:numpwc_T1D)])
    params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+(1:numpwc_T2D)])
    
    mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F]%*% params.reg_T2)
    
    res1 = obslogLn_pwc_vec_2(dt$R, dt$Y, dt$T0, dt$A1_dagger, dt$T2_dagger, dt$delta2, dt$T_dagger, dt$deltaD, 
                              dt$X1, brks_T2, params0_T2, beta1, mu_X2_T2, 
                              brks_T1D, params0_T1D, eta1, mu_X2_T1D, 
                              brks_T2D, params0_T2D, eta2, mu_X2_T2D, mu_X)
    res = -sum(res1)
    return(res)
  }
  
  dim <- (numpwc_T2+numpwc_T1D+numpwc_T2D) + length(main.cov)+1 #+length(TD.cov)
  est.res <- nlm(f=obs.f, p=c(rep(log(0.01), numpwc_T2), rep(0.0, length(main.cov)), 
                              rep(log(0.01), numpwc_T1D), rep(log(0.01), numpwc_T2D), 0.0),
                 hessian=T)
  
  vartheta = est.res$est
  params0_T2 = exp(vartheta[1:numpwc_T2])
  beta1 = 0
  params.reg_T2 = c(vartheta[numpwc_T2+(1:length(main.cov))])
  params0_T1D = exp(vartheta[numpwc_T2+length(main.cov)+(1:numpwc_T1D)])
  params0_T2D = exp(vartheta[numpwc_T2+length(main.cov)+numpwc_T1D+(1:numpwc_T2D)])
  
  mu_X2_T2 = as.vector(Xmat[,main.cov, drop=F] %*% params.reg_T2)
  
  logsurvF = survF1_g_A1_X2_pwc_vec(dt$T2_dagger, dt$A1_dagger, brks_T2, params0_T2, mu_X2_T2, 1)

  res = dt$delta2 + logsurvF

  return(res)
}


grs_noX1_noX2.f <- function(n2_samp, strata.n, dt, ref.group, score_res){
  
  score_res0 = score_res
  score_res = score_res0 - mean(score_res)
  
  order_resi = order(score_res)
  n = sum(strata.n)
  s_residual <- n2_samp*0.3
  s_id = (1:n)[c(order_resi[1:s_residual], order_resi[(n-(n2_samp-s_residual)+1):n])]
  mart.opt = score_res[s_id]
  best_var = var(mart.opt)
  for (k in (s_residual+1):(n2_samp*0.7)){
    s_id = (1:n)[c(order_resi[1:k], order_resi[(n-(n2_samp-k)+1):n])]
    mart.opt = score_res[s_id]
    tmp_var = var(mart.opt)
    
    #print(tmp_var)
    
    if (tmp_var > best_var) {
      s_residual = k
      best_var = tmp_var
    }
    
  }
  
  s_id = (1:n)[c(order_resi[1:s_residual], order_resi[(n-(n2_samp-s_residual)+1):n])]
  
  s_id = dt$id[s_id]
  
  s_dt = subset(dt, id %in% s_id)
  s_prob <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]/strata.n[g]
      
    }else{
      
      0
    }
    
  })
  
  res <- list(s_prob=s_prob,
              s_id=s_id)
  
  return(res)
  
}

grs_noX1.f <- function(n2_smp, strata.n, dt, ref.group, p_1_g_0, p_1_g_1, score_res){
  
  score_res0 = score_res
  indX2_1 <- which(dt$X2==1); indX2_0 <- which(dt$X2==0)
  
  score_res = score_res0 - mean(score_res0)
  
  score_res[indX2_0] = (score_res[indX2_0] )*sqrt(p_1_g_0*(1-p_1_g_0))
  score_res[indX2_1] = (score_res[indX2_1] )*sqrt(p_1_g_1*(1-p_1_g_1))
  order_resi0 = order(score_res[indX2_0])
  order_resi1 = order(score_res[indX2_1])
  
  nX2_0 <- length(indX2_0); nX2_1 <- length(indX2_1)
  
  s_residual = min(30, n2_smp/6, nX2_1/6)#min(n2_smp, nX2_1)/10#/4

  s_id = c(indX2_1[order_resi1[1:s_residual]], indX2_1[order_resi1[(nX2_1-s_residual+1):nX2_1]], 
           indX2_0[order_resi0[1:(n2_smp/2-s_residual)]], indX2_0[order_resi0[(nX2_0-(n2_smp/2-s_residual)+1):nX2_0]])

  mart.opt = score_res[s_id]
  X2.opt = dt$X2[s_id]
  Y.opt = dt$Y[s_id]
  
#  print(s_residual)
#  print(length(unique(s_id)))
#  print(head(s_id))
#  print(mart.opt)

  best_var = (var(mart.opt[which(X2.opt == 0)])*sum(X2.opt==0)+var(mart.opt[which(X2.opt == 1)])*sum(X2.opt==1))/n2_smp
  
  #print(head(which(X2.opt == 0)))  
  #print(s_residual)
  for (k in (s_residual+1):((min(n2_smp, nX2_1))/2-min(60, n2_smp/3, nX2_1/3))){
    
    s_id = c(indX2_1[order_resi1[1:k]], indX2_1[order_resi1[(nX2_1-k+1):nX2_1]], 
             indX2_0[order_resi0[1:(n2_smp/2-k)]], 
             indX2_0[order_resi0[(nX2_0-(n2_smp/2-k)+1):nX2_0]])
    
    mart.opt = score_res[s_id]
    X2.opt = dt$X2[s_id]
#    print(best_var)    
    tmp_var = (var(mart.opt[which(X2.opt == 0)])*sum(X2.opt==0)+var(mart.opt[which(X2.opt == 1)])*sum(X2.opt==1))/n2_smp
    
    if (tmp_var > best_var) {
      s_residual = k
      best_var = tmp_var
    }
    
  }
  
  #print(s_residual)
  s_id = c(indX2_1[order_resi1[1:s_residual]], indX2_1[order_resi1[(nX2_1-s_residual+1):nX2_1]], 
           indX2_0[order_resi0[1:(n2_smp/2-s_residual)]], indX2_0[order_resi0[(nX2_0-(n2_smp/2-s_residual)+1):nX2_0]])
  
  s_id = dt$id[s_id]
  
  s_dt = subset(dt, id %in% s_id)
  s_prob <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]/strata.n[g]
      
    }else{
      
      0
    }
    
  })
  
  s_m <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]
      
    }else{
      
      0
    }
    
  })
  
  res <- list(s_prob=s_prob,
              s_id=s_id,
              s_m=s_m)
  
  return(res)
  
}

rs_noX1.f <- function(n2_smp, strata.n, dt, ref.group, score_res){
  
  score_res0 = score_res
  indX2_1 <- which(dt$X2==1); indX2_0 <- which(dt$X2==0)
  
  score_res = score_res0 - mean(score_res0)
  
  score_res[indX2_0] = (score_res[indX2_0] )#*sqrt(p_1_g_0*(1-p_1_g_0))
  score_res[indX2_1] = (score_res[indX2_1] )#*sqrt(p_1_g_1*(1-p_1_g_1))
  order_resi0 = order(score_res[indX2_0])
  order_resi1 = order(score_res[indX2_1])
  
  nX2_0 <- length(indX2_0); nX2_1 <- length(indX2_1)
  
  s_residual = 30#min(n2_smp, nX2_1)/10#/4
  s_id = c(indX2_1[order_resi1[1:s_residual]], indX2_1[order_resi1[(nX2_1-s_residual+1):nX2_1]], 
           indX2_0[order_resi0[1:(n2_smp/2-s_residual)]], indX2_0[order_resi0[(nX2_0-(n2_smp/2-s_residual)+1):nX2_0]])
  
  mart.opt = score_res[s_id]
  X2.opt = dt$X2[s_id]
  Y.opt = dt$Y[s_id]
  
  best_var = (var(mart.opt[which(X2.opt == 0)])*sum(X2.opt==0)+var(mart.opt[which(X2.opt == 1)])*sum(X2.opt==1))/n2_smp
  
  #print(s_residual)
  for (k in (s_residual+1):((min(n2_smp, nX2_1))/2-60)){
    
    s_id = c(indX2_1[order_resi1[1:k]], indX2_1[order_resi1[(nX2_1-k+1):nX2_1]], 
             indX2_0[order_resi0[1:(n2_smp/2-k)]], 
             indX2_0[order_resi0[(nX2_0-(n2_smp/2-k)+1):nX2_0]])
    
    mart.opt = score_res[s_id]
    X2.opt = dt$X2[s_id]
    
    tmp_var = (var(mart.opt[which(X2.opt == 0)])*sum(X2.opt==0)+var(mart.opt[which(X2.opt == 1)])*sum(X2.opt==1))/n2_smp
    
    if (tmp_var > best_var) {
      s_residual = k
      best_var = tmp_var
    }
    
  }
  
  #print(s_residual)
  s_id = c(indX2_1[order_resi1[1:s_residual]], indX2_1[order_resi1[(nX2_1-s_residual+1):nX2_1]], 
           indX2_0[order_resi0[1:(n2_smp/2-s_residual)]], indX2_0[order_resi0[(nX2_0-(n2_smp/2-s_residual)+1):nX2_0]])
  
  s_id = dt$id[s_id]
  
  s_dt = subset(dt, id %in% s_id)
  s_prob <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]/strata.n[g]
      
    }else{
      
      0
    }
    
  })
  
  s_m <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]
      
    }else{
      
      0
    }
    
  })
  
  res <- list(s_prob=s_prob,
              s_id=s_id,
              s_m=s_m)
  
  return(res)
  
}

rs_noX1_multiX2.f <- function(n_sample, strata.n, dt, ref.group, resd_res){
  
  score_res = resd_res
  
  order_resi = order(score_res)
  #N1 =  sum(score_res>=0)
  k1 = round(n_sample/2)
  nX2 <- length(dt$id)
  
  s_res_end = k1#min(k1, N1)
  s_res_start = n_sample - s_res_end
  
  s_id =  c(dt$id[order_resi[(nX2-s_res_end+1):nX2]], 
            dt$id[order_resi[1:s_res_start]])
  
  s_dt = subset(dt, id %in% s_id)
  s_prob <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]/strata.n[g]
      
    }else{
      
      0
    }
    
  })
  
  res <- list(s_prob=s_prob,
              s_id=s_id)
  
  return(res)
  
}


approx_Delta_beta1.f <- function(gradstep, dt_use, Xmat, main.cov, brks_vec, est.res){
  
  #hess = -est.res$hessian
  score_mat = ipw_weighted_Pscore.pwc.f(gradstep, dt_use, Xmat, main.cov, est.res$x, brks_vec)
  hess = -ipw_weighted_Phess.pwc.f(gradstep, dt_use, Xmat, main.cov, est.res$x, brks_vec)
    #ipw_weighted_score_noX2.mortality.f(gradstep, dt_use, est.res$par)#ipw_weighted_score.f(gradstep, dt_use, brks_vec, est.res$par) 
  dUdw_vec = sapply(1:(dim(score_mat)[1]), function(j){
    colSums(score_mat[-j,])
  })
  
  res = t(ginv(hess) %*% dUdw_vec)
  
  if(is.infinite(brks_vec[1])){
    pos = 3
  }else{
    pos = length(brks_vec) + 2
  }
  
  return(res[,pos])
}

# from https://github.com/T0ngChen/multiwave/blob/master/nwts/nwtsparany.r 
integer.neyman.w2 = function(n.strata, NS, sample.size, lower, upper){
  nc = max(upper + 1)
  s = 1:nc
  arr = matrix(rep(NS, each = nc)/sqrt(s*(s+1)), nrow = n.strata, byrow = T)
  arr.list = as.list(as.data.frame(t(arr)))
  for(i in 1:length(arr.list)){
    arr.list[[i]][c(1:(lower[i]), (upper[i]):nc)] = 0
  }
  arr = do.call(cbind, arr.list)
  rk <- order(arr, na.last=TRUE, decreasing=TRUE)[1:(sample.size - sum(lower))]
  re.om.zero = table(rk%/%nc + 1) 
  re.zero = lower #rep(0, n.strata)
  re.zero[1:n.strata %in% names(re.om.zero)] = re.zero[1:n.strata %in% names(re.om.zero)] + re.om.zero
  n_str = re.zero
  
  n_str
}

Neyman_alloc_influence.f <- function(inf_vec, n_sample, lb.str.n, ub.str.n, phaseI_res){
  
  dt_c = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  #n <- max(dt_c$id)
  
  sd_vec = numeric()
  for(i in 1:length(ref.group)){
    sd_vec = c(sd_vec, sd(inf_vec[ dt_c$group %in% (ref.group[i])]))
  }
  
  s_ipw = integer.neyman.w2(n.strata = length(strata.n), NS = sd_vec*strata.n, sample.size = n_sample, 
                            lower = lb.str.n, upper = ub.str.n)
  
  id_s <- unlist(sapply(which(s_ipw!=0), function(x){
    as.numeric(sample(as.character(dt_c$id[dt_c$group %in% ref.group[x]]),s_ipw[x]))
  }))
  
  #if(length(id_s)<n_sample){
  #  id_s <- c(id_s, as.numeric(sample(as.character(dt_c$id[!(dt_c$id %in% id_s)]), n_sample-length(id_s))))
  #}
  
  dt_use = dt_c[(dt_c$id %in% id_s),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  prob_s = s_num/strata.n
  
  res <- list(s_id=id_s,
              s_prob=prob_s,
              s_ipw=s_ipw)
  
  return(res)
}

Neyman_adpt_2wave.f <- function(brks_vec, n2_samp_a, n2_samp, pcuts, design, design.factor, 
                                phaseI_res, Xmat, main.cov){
  
  phIIa.res <- subopt.f(n2_samp_a, "bal", phaseI_res)
  #ipw_subopt_AY.f(n2_samp_a, design, phaseI_res)
  
  dt_a_use <- phaseI_res$dt_ext #dt_c
  dt_a_use[!(dt_a_use$id %in% phIIa.res$s_id), "X1"] <- NA
  dt_a_use[!(dt_a_use$id %in% phIIa.res$s_id), "R"] <- 0
  dt_a_use$probs <- sapply(dt_a_use$group, function(x){
    phIIa.res$s_prob[phaseI_res$ref.strata %in% x]
  })
  #est_a <- ipw_est.idm.f(subset(dt_a_use, R==1))
  est_a <- ipw_estPL.nleqslv.pwc.f(subset(dt_a_use, R==1), Xmat[dt_a_use$R==1,], main.cov, brks_vec)
  
  del_beta1_a <- approx_Delta_beta1.f(1e-06, subset(dt_a_use, R==1), Xmat[dt_a_use$R==1,], main.cov, brks_vec, est_a)
  
  str.n.a = phIIa.res$s_m
  
  if(any(str.n.a<=2)){
    warning(paste0("untrusted sd. of influential fct in", which(str.n.a<=2)))
  }
  
  inf_vec <- del_beta1_a#/(subset(dt_a_use, R==1)$probs)
  
  sd_vec = numeric()
  for(i in 1:length(phaseI_res$ref.strata)){
    sd_vec = c(sd_vec, sd(inf_vec[ subset(dt_a_use, R==1)$group %in% (phaseI_res$ref.strata[i])]))
  }
  
  s_ipw.b = integer.neyman.w2(n.strata = length(str.n.a), NS = sd_vec*(phaseI_res$strata.n), #(phaseI_res$strata.n-str.n.a),#
                              sample.size = n2_samp-n2_samp_a, 
                              lower = rep(0.0, length(str.n.a))*(phaseI_res$strata.n-str.n.a), 
                              upper = phaseI_res$strata.n-str.n.a)
  #s_ipw = integer.neyman.w2(n.strata = length(str.n.a), NS = sd_vec*(phaseI_res$strata.n),
  #                          sample.size = n2_samp, 
  #                          lower = pmax(rep(0.05, length(str.n.a))*phaseI_res$strata.n, str.n.a), 
  #                          upper = phaseI_res$strata.n)
  #s_ipw.b=s_ipw - str.n.a
  dt_b=subset(dt_a_use, R==0)
  s_id.b <- unlist(sapply(which(s_ipw.b!=0), function(x){
    as.numeric(sample(as.character(dt_b$id[dt_b$group %in% phaseI_res$ref.strata[x]]),s_ipw.b[x]))
  }))
  
  #if(length(s_id.b)<(n2_samp-n2_samp_a)){
  #  s_id.b <- c(s_id.b, as.numeric(sample(as.character(dt_b$id[!(dt_b$id %in% s_id.b)]), (n2_samp-n2_samp_a)-length(s_id.b))))
  #}
  
  s_id = c(phIIa.res$s_id, s_id.b)
  dt_use = subset(phaseI_res$dt_ext, id %in% s_id) #dt_ext[dt_ext$id %in% s_id,]
  s_num <- sapply(phaseI_res$ref.strata, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  
  #print(sum(s_num))
  s.prob <- s_num/phaseI_res$strata.n
  #print(str.n.a)
  
  res = list(s_id=s_id,
             s_id.a=phIIa.res$s_id,
             s_id.b=s_id.b,
             s_prob = s.prob, 
             s_m=str.n.a+s_ipw.b,
             s.a=str.n.a,
             s.b=s_ipw.b)
  
  return(res)
}


