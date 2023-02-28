

approx_Delta_beta1_dm.f <- function(gradstep, dt_use, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, vartheta, simplified){
  
  #hess = -est.res$hessian
  score_mat = ipw_weighted_score.pwc_dm.f(gradstep, vartheta, dt_use, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified)
    #ipw_weighted_Pscore.pwc.f(gradstep, dt_use, Xmat, main.cov, est.res$x, brks_vec)
  hess = -ipw_weighted_hess.pwc_dm.f(gradstep, vartheta, dt_use, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified)
    #-ipw_weighted_Phess.pwc.f(gradstep, dt_use, Xmat, main.cov, est.res$x, brks_vec)
  #ipw_weighted_score_noX2.mortality.f(gradstep, dt_use, est.res$par)#ipw_weighted_score.f(gradstep, dt_use, brks_vec, est.res$par) 
  dUdw_vec = sapply(1:(dim(score_mat)[1]), function(j){
    colSums(score_mat[-j,])
  })
  
  res = t(ginv(hess) %*% dUdw_vec)
  
  if(is.infinite(brks_T2[1])){
    pos = 2
  }else{
    pos = length(brks_T2) + 2
  }
  
  return(res[,pos])
}


Neyman_adpt_2wave_dm.f <- function(gradstep, brks_T2, brks_T1D, brks_T2D, n2_samp_a, n2_samp, 
                                   pcuts, phI.design.factor, #phII.design, 
                                   phaseI_res, Xmat, main.cov, simplified){
  
  phIIa.res <- subopt.f(n2_samp_a, "bal", phaseI_res)
  #ipw_subopt_AY.f(n2_samp_a, design, phaseI_res)
  
  dt_a_use <- phaseI_res$dt_ext #dt_c
  dt_a_use[!(dt_a_use$id %in% phIIa.res$s_id), "X1"] <- NA
  dt_a_use[!(dt_a_use$id %in% phIIa.res$s_id), "R"] <- 0
  dt_a_use$probs <- sapply(dt_a_use$group, function(x){
    phIIa.res$s_prob[phaseI_res$ref.strata %in% x]
  })
  #est_a <- ipw_est.idm.f(subset(dt_a_use, R==1))
  #est_a <- ipw_estPL.nleqslv.pwc.f(subset(dt_a_use, R==1), Xmat[dt_a_use$R==1,], main.cov, brks_vec)
  est_a <- ipw_estL.pwc_dm.f(subset(dt_a_use, R==1), Xmat[dt_a_use$R==1,,drop=F], main.cov, brks_T2, brks_T1D, brks_T2D, simplified)
  
  del_beta1_a <- approx_Delta_beta1_dm.f(gradstep, subset(dt_a_use, R==1), Xmat[dt_a_use$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, est_a$par, simplified)
  
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
  
  s_id = c(phIIa.res$s_id, s_id.b)
  dt_use = subset(phaseI_res$dt_ext, id %in% s_id) #dt_ext[dt_ext$id %in% s_id,]
  s_num <- sapply(phaseI_res$ref.strata, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  
  print(sum(s_num))
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





