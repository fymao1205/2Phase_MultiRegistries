

rldt_ml_script_DM.f <- function(nsim, n2_samp, prop.a=NULL, dt, Xmat, design, design.factor=c("Y", "del2"), pcuts,
                             main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
  #pcuts = c(0.33,0.66)
  phaseI.res <- phase_I_stratify_IDM.f(pcuts, dt, design.factor)
  
  phI.strat.style = paste0(design.factor, collapse = "_")
  if(simplified){
    file1.out <- paste(paste0(n2_samp, "ML.", phI.strat.style, paste0("pcuts", pcuts, collapse=""), design, prop.a, 
                              paste0("main.X",main.cov, collapse = ""), 
                              paste0("nui.X",nui.cov, collapse = ""), 
                              paste0("brks_T2", brks_T2, collapse = ""),
                              paste0("brks_T1D", brks_T1D, collapse = ""),
                              paste0("phi", collapse = ""),
                              ".dat", collapse=""), sep="")
  }else{
    file1.out <- paste(paste0(n2_samp, "ML.", phI.strat.style, paste0("pcuts", pcuts, collapse=""), design, prop.a, 
                              paste0("main.X",main.cov, collapse = ""), 
                              paste0("nui.X",nui.cov, collapse = ""), 
                              paste0("brks_T2", brks_T2, collapse = ""),
                              paste0("brks_T1D", brks_T1D, collapse = ""),
                              paste0("brks_T2D", brks_T2D, collapse = ""),
                              ".dat", collapse=""), sep="")
  }
  
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  
  if(any(is.finite(brks_T2))){
    num.pwc_T2 = length(brks_T2)+1
  }else{
    num.pwc_T2=1
  }
  
  if(any(is.finite(brks_T1D))){
    num.pwc_T1D = length(brks_T1D)+1
  }else{
    num.pwc_T1D=1
  }
  
  if(!simplified){
    if(any(is.finite(brks_T2D))){
      num.pwc_T2D = length(brks_T2D)+1
    }else{
      num.pwc_T2D=1
    }
  }else{
    num.pwc_T2D=1
  }
  
  num.eta <- ifelse(is.null(nui.cov), 1, 1+length(nui.cov))
  if(simplified){
    label.tx <- c("nsim", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                            paste0("logpsi1.", 1:num.pwc_T1D), 
                            paste0("phi"), 
                            paste0("eta", (1:num.eta)-1)),
                  paste0("ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                   paste0("logpsi1.", 1:num.pwc_T1D), 
                                   paste0("phi"), 
                                   paste0("eta", (1:num.eta)-1))), 
                  paste0("sw.ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                      paste0("logpsi1.", 1:num.pwc_T1D), 
                                      paste0("phi"), 
                                      paste0("eta", (1:num.eta)-1))),
                  paste0("num.score.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                         paste0("logpsi1.", 1:num.pwc_T1D), 
                                         paste0("phi"), 
                                         paste0("eta", (1:num.eta)-1))),
                  paste0("g", phaseI.res$ref.strata), paste0("sprob.g", phaseI.res$ref.strata),
                  paste0("sY", 1:2))
  }else{
    label.tx <- c("nsim", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                            paste0("logpsi1.", 1:num.pwc_T1D), 
                            paste0("logpsi2.", 1:num.pwc_T2D), 
                            paste0("eta", (1:num.eta)-1)),
                  paste0("ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                   paste0("logpsi1.", 1:num.pwc_T1D), 
                                   paste0("logpsi2.", 1:num.pwc_T2D), 
                                   paste0("eta", (1:num.eta)-1))), 
                  paste0("sw.ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                      paste0("logpsi1.", 1:num.pwc_T1D), 
                                      paste0("logpsi2.", 1:num.pwc_T2D), 
                                      paste0("eta", (1:num.eta)-1))),
                  paste0("num.score.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                         paste0("logpsi1.", 1:num.pwc_T1D), 
                                         paste0("logpsi2.", 1:num.pwc_T2D), 
                                         paste0("eta", (1:num.eta)-1))),
                  paste0("g", phaseI.res$ref.strata), paste0("sprob.g", phaseI.res$ref.strata),
                  paste0("sY", 1:2))
  }
  
  cat(as.character(label.tx), sep=" ", "\n", append=T, file=file1.out)
  
  iter=1
  for(m in 1:(2*nsim)){
    
    if(iter>nsim){break}
    
    set.seed(m)
    
    if(design=="tao"){
      
      Mmu <- coxph(Surv(T2_dagger, delta2)~a1.ind.g20.l40+a1.ind.g40+X2, phaseI.res$dt_ext)$resid
      if(isTRUE(class(Mmu)=="try-error")) { next }
      s1.res <- rs_noX1.f(n2_samp, strata.n=phaseI.res$strata.n, phaseI.res$dt_ext, ref.group=phaseI.res$ref.strata, score_res=Mmu)
      
    }else if(design=="ext_Mmu"){
      
      #Mmu <- scorePL_mu_pwc.f(phaseI.res$dt_ext, Xmat, main.cov, brks)
      Mmu <- try(scoreL_mu_pwc_dm.f(phaseI.res$dt_ext, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
      if(isTRUE(class(Mmu)=="try-error")) { next }
      s1.res <- rs_noX1_multiX2.f(n2_samp, strata.n=phaseI.res$strata.n, phaseI.res$dt_ext, ref.group=phaseI.res$ref.strata, Mmu)
      
    }else if(design=="awext_Mmu"){
      
      ## phase IIa
      IIa.res <- try(subopt.f(n2_samp*prop.a, design="srs", phaseI.res))
      dt_IIa <- phaseI.res$dt_ext 
      dt_IIa[!(dt_IIa$id %in% IIa.res$s_id), "X1"] <- NA
      dt_IIa[!(dt_IIa$id %in% IIa.res$s_id), "R"] <- 0
      dt_IIa$probs <- sapply(dt_IIa$group, function(x){
        IIa.res$s_prob[x == phaseI.res$ref.strata]
      })
      
      sIIa.est <- try(estL.pwc_dm.f(dt_IIa, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified))
      if(isTRUE(class(sIIa.est)=="try-error")) { 
        next
      }
      if(is.null(nui.cov)){
        mu_X = as.vector(rep(c(sIIa.est$est[num.pwc_T2+length(main.cov)+num.pwc_T1D+num.pwc_T2D+2]), dim(Xmat)[1]))
      }else{
        params.X = c(sIIa.est$est[num.pwc_T2+length(main.cov)+num.pwc_T1D+num.pwc_T2D+1+(1:(length(nui.cov)+1))])
        mu_X = as.vector(params.X[1]+Xmat[,nui.cov, drop=F] %*% params.X[-1])
      }
      p_1_g_x2 <- expit(mu_X)
      var_x1_g_x2 <- p_1_g_x2*(1-p_1_g_x2)
      
      phaseI_res.b = phaseI.res
      phaseI_res.b$strata.n <- phaseI.res$strata.n*(1-IIa.res$s_prob)
      phaseI_res.b$Abar_cut <- phaseI.res$Abar_cut
      phaseI_res.b$dt_ext <- subset(phaseI.res$dt_ext, !(id %in% IIa.res$s_id))
      phaseI_res.b$ref.strata <- phaseI.res$ref.strata
      nsmp.b=n2_samp*(1-prop.a)
      
      Mmu <- try(scoreL_mu_pwc_dm.f(phaseI.res$dt_ext, Xmat, main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
      if(isTRUE(class(Mmu)=="try-error")) { next }
      var_x1_g_x2.b <- var_x1_g_x2[!(dt_IIa$id %in% IIa.res$s_id)]
      Mmu.b <- Mmu[!(dt_IIa$id %in% IIa.res$s_id)]*sqrt(var_x1_g_x2.b)
      s0.res <- rs_noX1_multiX2.f(n2_samp*(1-prop.a), strata.n=phaseI_res.b$strata.n, phaseI_res.b$dt_ext, 
                           ref.group=phaseI_res.b$ref.strata, Mmu.b)
      s1.res <- s0.res
      s1.res$s_id <- c(IIa.res$s_id, s0.res$s_id)
      s1.res$s_m <- c(IIa.res$s_m+s0.res$s_m)
      s1.res$sprob <- s1.res$s_m/phaseI.res$strata.n
      
    }else{
      
      s1.res <- subopt.f(n2_samp, design, phaseI.res)
      
    }
    
    dt_s1 <- phaseI.res$dt_ext
    dt_s1[!(dt_s1$id %in% s1.res$s_id), "X1"] <- NA
    dt_s1[!(dt_s1$id %in% s1.res$s_id), "R"] <- 0
    #s1.est <- try(estPL.pwc.f(dt_s1, Xmat, main.cov, nui.cov, brks))
    s1.est <- try(estL.pwc_dm.f(dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified))
    if(isTRUE(class(s1.est)=="try-error")) { 
      next
    }
    #s1.ase = sw.asePL.pwc.f(1e-03, s1.est$par, dt_s1, Xmat, main.cov, nui.cov, brks) #sw.asePL.pwc.f(1e-06, pl.est.res$par, dt, Xmat, main.cov, nui.cov, brks_T2)
    s1.ase = try(sw.aseL.pwc_dm.f(1e-06, s1.est$est, dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D, simplified))
    
    s.Y = as.data.frame(xtabs(~Y, subset(dt_s1, R==1)))$Freq
    
    cat(m, s1.est$est, sqrt(diag(ginv(s1.est$hess))), 
        s1.ase, s1.res$s_prob, s1.est$gradient, 
        s1.res$s_m, s.Y,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
  
}


rldt_ipw_script_DM.f <- function(nsim, n2_samp, prop.a=NULL, dt, Xmat, design, design.factor=c("Y", "del2"), pcuts,
                                main.cov, brks_T2, brks_T1D, brks_T2D, simplified){
  
  #pcuts = c(0.33,0.66)
  phaseI.res <- phase_I_stratify_IDM.f(pcuts, dt, design.factor)
  
  phI.strat.style = paste0(design.factor, collapse = "_")
  if(simplified){
    file1.out <- paste(paste0(n2_samp, "IPW.", phI.strat.style, paste0("pcuts", pcuts, collapse=""), 
                              design, prop.a, paste0("X",main.cov, collapse = ""), 
                              paste0("brks_T2", brks_T2, collapse = ""),
                              paste0("brks_T1D", brks_T1D, collapse = ""),
                              paste0("phi"),
                              ".dat", collapse=""), sep="")
  }else{
    file1.out <- paste(paste0(n2_samp, "IPW.", phI.strat.style, paste0("pcuts", pcuts, collapse=""), 
                              design, prop.a, paste0("X",main.cov, collapse = ""), 
                              paste0("brks_T2", brks_T2, collapse = ""),
                              paste0("brks_T1D", brks_T1D, collapse = ""),
                              paste0("brks_T2D", brks_T2D, collapse = ""),
                              ".dat", collapse=""), sep="")
  }
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  
  if(any(is.finite(brks_T2))){
    num.pwc_T2 = length(brks_T2)+1
  }else{
    num.pwc_T2=1
  }
  
  if(any(is.finite(brks_T1D))){
    num.pwc_T1D = length(brks_T1D)+1
  }else{
    num.pwc_T1D=1
  }
  
  if(!simplified){
    if(any(is.finite(brks_T2D))){
      num.pwc_T2D = length(brks_T2D)+1
    }else{
      num.pwc_T2D=1
    }
  }
  
  #num.eta <- ifelse(is.null(nui.cov), 1, 1+length(nui.cov))
  if(simplified){
    label.tx <- c("nsim", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                            paste0("logpsi1.", 1:num.pwc_T1D), 
                            paste0("phi")),
                  paste0("ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                   paste0("logpsi1.", 1:num.pwc_T1D), 
                                   paste0("phi"))), 
                  paste0("num.score.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                         paste0("logpsi1.", 1:num.pwc_T1D), 
                                         paste0("phi"))),
                  paste0("g", phaseI.res$ref.strata), paste0("sprob.g", phaseI.res$ref.strata),
                  paste0("sY", 1:2))
  }else{
    label.tx <- c("nsim", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                            paste0("logpsi1.", 1:num.pwc_T1D), 
                            paste0("logpsi2.", 1:num.pwc_T2D)),
                  paste0("ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                   paste0("logpsi1.", 1:num.pwc_T1D), 
                                   paste0("logpsi2.", 1:num.pwc_T2D))), 
                  #paste0("sw.ASE.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                  #                    paste0("logpsi1.", 1:num.pwc_T1D), 
                  #                    paste0("logpsi2.", 1:num.pwc_T2D))),
                  paste0("num.score.", c(paste0("logalp.", 1:num.pwc_T2), paste0("beta", c(1, 1+main.cov)), 
                                         paste0("logpsi1.", 1:num.pwc_T1D), 
                                         paste0("logpsi2.", 1:num.pwc_T2D))),
                  paste0("g", phaseI.res$ref.strata), paste0("sprob.g", phaseI.res$ref.strata),
                  paste0("sY", 1:2))
  }
  
  cat(as.character(label.tx), sep=" ", "\n", append=T, file=file1.out)
  
  iter=1
  for(m in 1:(2*nsim)){
    
    if(iter>nsim){break}
    
    set.seed(m)
    
    if(design=="neyA"){
      
      n2_samp_a = n2_samp*prop.a
      s1.res <- Neyman_adpt_2wave_dm.f(gradstep=1e-06, brks_T2, brks_T1D, brks_T2D, n2_samp_a, n2_samp, 
                                       pcuts, phI.design.factor, #phII.design, 
                                       phaseI.res, Xmat, main.cov, simplified)
      
    }else{
      
      s1.res <- subopt.f(n2_samp, design, phaseI.res)
      
    }
    
    dt_s1 <- phaseI.res$dt_ext
    dt_s1[!(dt_s1$id %in% s1.res$s_id), "X1"] <- NA
    dt_s1[!(dt_s1$id %in% s1.res$s_id), "R"] <- 0
    dt_s1$probs <- sapply(dt_s1$group, function(x){
      s1.res$s_prob[x == phaseI.res$ref.strata]
    })
    #s1.est <- try(estL.pwc_dm.f(dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    s1.est <- try(ipw_estL.pwc_dm.f(subset(dt_s1, R==1), Xmat[dt_s1$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
    if(isTRUE(class(s1.est)=="try-error")) { 
      next
    }
    
    #s1.ase = try(sw.aseL.pwc_dm.f(1e-06, s1.est$est, dt_s1, Xmat, main.cov, nui.cov, brks_T2, brks_T1D, brks_T2D))
    
    s1.ase <- try(ipw_estL.pwc.ase_dm.f(1e-06, s1.est$par, phaseI.res$ref.strata, subset(dt_s1, R==1), Xmat[dt_s1$R==1,], main.cov, brks_T2, brks_T1D, brks_T2D, simplified))
    if(isTRUE(class(s1.ase)=="try-error")) { 
      next
    }
    
    s.Y = as.data.frame(xtabs(~Y, subset(dt_s1, R==1)))$Freq
    
    cat(m, s1.est$par, s1.ase$se, s1.ase$score,  
        s1.res$s_m, s1.res$s_prob, s.Y,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
  
}

