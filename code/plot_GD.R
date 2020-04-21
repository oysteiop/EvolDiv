


plot_GD=function(gmat=out$G, dmat=out$D, pmat=NULL, plot=FALSE, species="Species", xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL){

  res = data.frame(traits = c("Original", "Original","G eigenvectors", "D eigenvectors", "D eigenvectors", "P eigenvectors", "P eigenvectors", "All"),
                   measure = c("e","c","e","e","c","e","c","e"),
                   slope = NA,
                   r2 = NA)
  
  # Compute eigenvectors etc.
  g_ev = eigen(gmat)$vectors
  var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
  var_g_g_c = evolvabilityBeta(gmat, Beta = g_ev)$c
  var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
  
  d_ev = eigen(dmat)$vectors
  var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
  var_g_d_c = evolvabilityBeta(gmat, Beta = d_ev)$c
  var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e
  
  if(!is.null(pmat)){
  p_ev = eigen(pmat)$vectors
  var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
  var_g_p_c = evolvabilityBeta(gmat, Beta = p_ev)$c
  var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e
  }
  
  # Theta  
  theta = acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
  
  # Compute summary stats
  cvals = NULL
  for(i in 1:ncol(gmat)){
    b = rep(0, ncol(gmat))
    b[i] = 1
    cvals[i] = evolvabilityBeta(gmat, b)$c
  }
  
  mt = lm(log(diag(dmat))~log(diag(gmat)))
  beta_t = summary(mt)$coef[2,1]
  r2_t = summary(mt)$r.squared
  
  mtc = lm(log(diag(dmat))~log(cvals))
  beta_tc = summary(mtc)$coef[2,1]
  r2_tc = summary(mtc)$r.squared
  
  mg = lm(log(var_d_g)~log(var_g_g))
  beta_g = summary(mg)$coef[2,1]
  r2_g = summary(mg)$r.squared
  
  md = lm(log(var_d_d)~log(var_g_d))
  beta_d = summary(md)$coef[2,1]
  r2_d = summary(md)$r.squared
  
  mdc = lm(log(var_d_d)~log(var_g_d_c))
  beta_dc = summary(mdc)$coef[2,1]
  r2_dc = summary(mdc)$r.squared
  
  ma = lm(log(c(diag(dmat), var_d_g, var_d_d))~log(c(diag(gmat), var_g_g, var_g_d)))
  beta_a = summary(ma)$coef[2,1]
  r2_a = summary(ma)$r.squared
  
  if(!is.null(pmat)){
  mp = lm(log(var_d_p)~log(var_g_p))
  beta_p = summary(mp)$coef[2,1]
  r2_p = summary(mp)$r.squared
  
  mpc = lm(log(var_d_p)~log(var_g_p_c))
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
     main=species, las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
if(!is.null(pmat)){
  points(log10(var_g_p), log10(var_d_p), pch=16, col="firebrick")
  legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors", "P eigenvectors"), 
         pch=c(16, 1, 16, 16), col=c("blue3", "black", "black", "firebrick"))
}
if(is.null(pmat)){
  legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors"), 
         pch=c(16, 1, 16), col=c("blue3", "black", "black"))
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
       main=species, las=1)
  points(log10(var_g_d_c), log10(var_d_d), pch=16)
  
  points(log10(cvals), log10(diag(dmat)), pch=16, col="blue3")
  if(!is.null(pmat)){
    points(log10(var_g_p_c), log10(var_d_p), pch=16, col="firebrick")
    legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors", "P eigenvectors"), 
           pch=c(16, 1, 16, 16), col=c("blue3", "black", "black", "firebrick"))
  }
  if(is.null(pmat)){
    legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors"), 
           pch=c(16, 1, 16), col=c("blue3", "black", "black"))
  }
}

outlist = list(res=res, theta = theta)

return(outlist)
}
