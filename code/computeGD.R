check=F
if(check){
  gmat=gmat
  dmat=dmat
  pmat=MeanP
  plot=FALSE 
  species="Species"
  nFam=out$nFam 
  nPop=out$nPop 
  SE=TRUE
  nSample=10 
  fixD=TRUE
  xmin=NULL 
  xmax=NULL 
  ymin=NULL 
  ymax=NULL
}


computeGD=function(gmat=out$G, dmat=out$D, pmat=NULL, plot=FALSE, species="Species",
                   nFam=out$nFam, nPop=out$nPop, SE=FALSE, nSample=100, fixD=TRUE,
                   xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL){
  
  res = data.frame(traits = c("Original", "Original","G eigenvectors", "D eigenvectors", "D eigenvectors", "P eigenvectors", "P eigenvectors", "All"),
                   measure = c("e","c","e","e","c","e","c","e"),
                   slope = NA,
                   slope_MC = NA,
                   SE = NA,
                   r2 = NA)
  
  dpost = dmat
  
  if(!isSymmetric.matrix(dmat)){
    dmat = matrix(apply(dpost, 2, median), nrow=sqrt(ncol(as.matrix(dpost))))
  }
  
  a = alignMat(gmat, dmat)
  
  # Compute eigenvectors etc.
  g_ev = eigen(a$G_ge)$vectors
  var_g_g = evolvabilityBeta(a$G_ge, Beta = g_ev)$e
  var_g_g_c = evolvabilityBeta(a$G_ge, Beta = g_ev)$c
  var_d_g = evolvabilityBeta(a$D_ge, Beta = g_ev)$e
  
  d_ev = eigen(a$D_de)$vectors
  var_g_d = evolvabilityBeta(a$G_de, Beta = d_ev)$e
  var_g_d_c = evolvabilityBeta(a$G_de, Beta = d_ev)$c
  var_d_d = evolvabilityBeta(a$D_de, Beta = d_ev)$e
  
  if(!is.null(pmat)){
    p_ev = eigen(pmat)$vectors
    var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
    var_g_p_c = evolvabilityBeta(gmat, Beta = p_ev)$c
    var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e
  }
  
  # Theta  
  theta = acos(t(eigen(gmat)$vectors[,1]) %*% eigen(dmat)$vectors[,1])*(180/pi)
  
  # Compute summary stats
  cvals = NULL
  for(i in 1:ncol(gmat)){
    b = rep(0, ncol(gmat))
    b[i] = 1
    cvals[i] = evolvabilityBeta(gmat, b)$c
    #cvals[i] = evolvabilityBeta(bendMat(gmat), b)$c
  }
  
  mt = lm(log(diag(dmat))~log(diag(gmat)), na=na.exclude)
  beta_t = summary(mt)$coef[2,1]
  r2_t = summary(mt)$r.squared
  
  mtc = lm(log(diag(dmat))~log(cvals), na=na.exclude)
  beta_tc = summary(mtc)$coef[2,1]
  r2_tc = summary(mtc)$r.squared
  
  mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
  beta_g = summary(mg)$coef[2,1]
  r2_g = summary(mg)$r.squared
  
  md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
  beta_d = summary(md)$coef[2,1]
  r2_d = summary(md)$r.squared
  
  mdc = lm(log(var_d_d)~log(var_g_d_c), na=na.exclude)
  beta_dc = summary(mdc)$coef[2,1]
  r2_dc = summary(mdc)$r.squared
  
  ma = lm(log(c(diag(dmat), var_d_g, var_d_d))~log(c(diag(gmat), var_g_g, var_g_d)), na=na.exclude)
  beta_a = summary(ma)$coef[2,1]
  r2_a = summary(ma)$r.squared
  
  mac = lm(log(c(diag(dmat), var_d_g, var_d_d))~log(c(cvals, var_g_g_c, var_g_d_c)), na=na.exclude)
  beta_ac = summary(mac)$coef[2,1]
  r2_ac = summary(mac)$r.squared
  
  if(!is.null(pmat)){
    mp = lm(log(var_d_p)~log(var_g_p), na=na.exclude)
    beta_p = summary(mp)$coef[2,1]
    r2_p = summary(mp)$r.squared
    
    mpc = lm(log(var_d_p)~log(var_g_p_c), na=na.exclude)
    beta_pc = summary(mpc)$coef[2,1]
    r2_pc = summary(mpc)$r.squared
    
    #ma = lm(log(c(diag(dmat), var_d_g, var_d_d, var_d_p))~log(c(diag(gmat), var_g_g, var_g_d, var_g_p)))
    #beta_a = summary(ma)$coef[2,1]
    #r2_a = summary(ma)$r.squared
  }
  
  if(is.null(pmat)){
    res$slope[c(1:5, 8)] = round(c(beta_t, beta_tc, beta_g, beta_d, beta_dc, beta_a), 3)
    res$r2[c(1:5, 8)] = round(c(r2_t, r2_tc, r2_g, r2_d, r2_dc, r2_a), 3)
  }
  if(!is.null(pmat)){
    res$slope = round(c(beta_t, beta_tc, beta_g, beta_d, beta_dc, beta_p, beta_pc, beta_a), 3)
    res$r2 = round(c(r2_t, r2_tc, r2_g, r2_d, r2_dc, r2_p, r2_pc, r2_a), 3)
  }
  
  # Generating standard errors by resampling
  if(SE){
    print("Computing standard errors by resampling")
    require(mvtnorm)
    source("code/bendMat.R")
    sbeta_g = sbeta_t = sbeta_tc = sbeta_d = sbeta_dc = sbeta_p = sbeta_pc = sbeta_a = NULL
    for(i in 1:nSample){
      if(i==nSample/2){print("Halfways done!")}
      sgmat = round(cov(rmvnorm(nFam, mean=rep(1, ncol(gmat)), sigma=gmat)), 10)
      #sdmat = round(cov(rmvnorm(nPop, mean=rep(1, ncol(dmat)), sigma=dmat)), 10)
      sdmat = matrix(dpost[sample(1:nSample, 1),], nrow=ncol(dmat))
      if(fixD){sdmat = dmat} 
      
      sa = alignMat(sgmat, sdmat)
      
      smt = lm(log(diag(sdmat))~log(diag(sgmat)), na=na.exclude)
      sbeta_t[i] = summary(smt)$coef[2,1]
      
      scvals = NULL
      for(ci in 1:ncol(sgmat)){
        b = rep(0, ncol(sgmat))
        b[ci] = 1
        scvals[ci] = evolvabilityBeta(sgmat, b)$c
      }
      
      if(sum(scvals>0)>0){
      smtc = lm(log(diag(sdmat))~log(scvals), na=na.exclude)
      sbeta_tc[i] = summary(smtc)$coef[2,1]
      }
      
      sg_ev = eigen(sa$G_ge)$vectors
      svar_g_g = evolvabilityBeta(sa$G_ge, Beta = sg_ev)$e
      svar_g_g_c = evolvabilityBeta(sa$G_ge, Beta = sg_ev)$c
      svar_d_g = evolvabilityBeta(sa$D_ge, Beta = sg_ev)$e
      
      smg = lm(log(svar_d_g)~log(svar_g_g), na=na.exclude)
      sbeta_g[i] = summary(smg)$coef[2,1]
      
      sd_ev = eigen(sa$D_de)$vectors
      svar_g_d = evolvabilityBeta(sa$G_de, Beta = sd_ev)$e
      svar_g_d_c = evolvabilityBeta(sa$G_de, Beta = sd_ev)$c
      svar_d_d = evolvabilityBeta(sa$D_de, Beta = sd_ev)$e
      
      smd = lm(log(svar_d_d)~log(svar_g_d), na=na.exclude)
      sbeta_d[i] = summary(smd)$coef[2,1]
      
      smdc = lm(log(svar_d_d)~log(svar_g_d_c), na=na.exclude)
      sbeta_dc[i] = summary(smdc)$coef[2,1]
      
      sma = lm(log(c(diag(sdmat), svar_d_g, svar_d_d))~log(c(diag(sgmat), svar_g_g, svar_g_d)), na=na.exclude)
      sbeta_a[i] = summary(sma)$coef[2,1]
      
      if(!is.null(pmat)){
        p_ev = eigen(pmat)$vectors
        svar_g_p = evolvabilityBeta(sgmat, Beta = p_ev)$e
        svar_g_p_c = evolvabilityBeta(sgmat, Beta = p_ev)$c
        svar_d_p = evolvabilityBeta(sdmat, Beta = p_ev)$e
        
        smp = lm(log(svar_d_p)~log(svar_g_p), na=na.exclude)
        sbeta_p[i] = summary(smp)$coef[2,1]
        
        smpc = lm(log(svar_d_p)~log(svar_g_p_c), na=na.exclude)
        sbeta_pc[i] = summary(smpc)$coef[2,1]
      }
    }
    if(is.null(pmat)){
      res$slope_MC[c(1:5, 8)] = round(c(mean(sbeta_t), mean(sbeta_tc), mean(sbeta_g), mean(sbeta_d), mean(sbeta_dc), mean(sbeta_a)), 3)
      res$SE[c(1:5, 8)] = round(c(sd(sbeta_t), sd(sbeta_tc), sd(sbeta_g), sd(sbeta_d), sd(sbeta_dc), sd(sbeta_a)), 3)
    }
    if(!is.null(pmat)){
      res$slope_MC[c(1:8)] = round(c(mean(sbeta_t), mean(sbeta_tc), mean(sbeta_g), mean(sbeta_d), mean(sbeta_dc), mean(sbeta_p), mean(sbeta_pc), mean(sbeta_a)), 3)
      res$SE[1:8] = round(c(sd(sbeta_t), sd(sbeta_tc), sd(sbeta_g), sd(sbeta_d), sd(sbeta_dc), sd(sbeta_p), sd(sbeta_pc), sd(sbeta_a)), 3)
    }
  }
  
  # Plot
  if(plot=="e"){
    x11(width=5, height=5)
    
    if(is.null(xmin)){xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))}
    if(is.null(xmax)){xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))}
    if(is.null(ymin)){ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))}
    if(is.null(ymax)){ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))}
    plot(log10(var_g_g), log10(var_d_g), 
         xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
         xlab="log10 (Evolvability [%])", 
         ylab="log10 (Divergence [x100])", 
         pch=16, main=species, las=1)
    points(log10(var_g_d), log10(var_d_d), pch=1)
    points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue3")
    #lines(-100:100, lty=2)
    
    #legend("topleft", legend=substitute(paste(r^2, " = ", v), list(v=signif(r2_a, 2))), bty="n")
    
    mean1 = mean(log10(c(diag(dmat), var_d_g, var_d_d)))
    mean2 = mean(log10(c(diag(gmat), var_g_g, var_g_d)))
    segments(x0=mean2-10, y0=mean1-10, x1=mean2+10, y1=mean1+10)
    
    if(!is.null(pmat)){
      points(log10(var_g_p), log10(var_d_p), pch=16, col="firebrick")
      #legend("bottomright", c("Original traits", "G directions", "D directions", "P directions"), 
      #       pch=c(16, 1, 16, 16), col=c("blue3", "black", "black", "firebrick"))
      legend("bottomright", legend=c(paste("Original traits (", round(100*r2_t, 1),")"),
                                     paste("G directions (", round(100*r2_g, 1),")"),
                                     paste("D directions (", round(100*r2_d, 1),")"),
                                     paste("P directions (", round(100*r2_p, 1),")")),
             pch=c(16, 16, 1, 16), col=c("blue3", "black", "black", "firebrick"))
    }
    if(is.null(pmat)){
      legend("bottomright", c("Original traits", "G directions", "D directions"), 
             pch=c(16, 16, 1), col=c("blue3", "black", "black"))
    }
  }
  
  if(plot=="c"){
    x11(width=5, height=5)
    
    if(is.null(xmin)){xmin = log10(min(c(var_g_g_c, var_g_d_c), na.rm=T))}
    if(is.null(xmax)){xmax = log10(max(c(var_g_g_c, var_g_d_c), na.rm=T))}
    if(is.null(ymin)){ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))}
    if(is.null(ymax)){ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))}
    plot(log10(var_g_g_c), log10(var_d_g), 
         xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
         xlab="log10 (Conditional evolvability [%])", 
         ylab="log10 (Divergence [x100])", 
         pch=16, main=species, las=1)
    points(log10(var_g_d_c), log10(var_d_d), pch=1)
    
    points(log10(cvals), log10(diag(dmat)), pch=16, col="blue3")
    
    #legend("topleft", legend=substitute(paste(r^2, " = ", v), list(v=signif(r2_ac, 2))), bty="n")
    
    mean1 = mean(log10(c(diag(dmat), var_d_g, var_d_d)))
    mean2 = mean(log10(c(cvals, var_g_g_c, var_g_d_c)))
    segments(x0=mean2-10, y0=mean1-10, x1=mean2+10, y1=mean1+10)
    
    if(!is.null(pmat)){
      points(log10(var_g_p_c), log10(var_d_p), pch=16, col="firebrick")
      legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors", "P eigenvectors"), 
             pch=c(16, 16, 1, 16), col=c("blue3", "black", "black", "firebrick"))
    }
    if(is.null(pmat)){
      legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors"), 
             pch=c(16, 16, 1), col=c("blue3", "black", "black"))
    }
  }
  
  outlist = list(res=res, theta = theta, 
                 evolvabilityMeans = evolvabilityMeans(a$G_ge), 
                 divergenceMeans = evolvabilityMeans(a$D_de))
  
  return(outlist)
}
