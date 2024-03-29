computeDelta = function(G, means, z0, ndrift=1000){
  require(mvtnorm)
  outdat = matrix(NA, nrow=nrow(means), ncol=10)

  #Drift expectation
  rb = randomDelta(n = ndrift, gmat = G)
  e_beta = NULL
  e_drift = mean(evolvabilityBeta(G*100, rb)$e)

for(i in 1:nrow(means)){
  z1 = unlist(means[i,])
  delta = log(z1)-log(z0)
  scale_delta = delta/sqrt(sum(delta^2)) 
  d = delta
  div = mean(abs(d))*100
  
  e_delta = evolvabilityBeta(G*100, scale_delta)$e
  c_delta = evolvabilityBeta(G*100, scale_delta)$c
  theta = acos(t(eigen(G)$vectors[,1]) %*% scale_delta)*(180/pi)
  
  outdat[i,]=c(div, e_delta, c_delta, theta, evolvabilityMeans(G*100)[c(1:4,7)], e_drift)
  }
return(outdat)
}

randomDelta = function(n = 1, gmat = NULL){
  X = t(rmvnorm(n, sigma=gmat))
  X = t(t(X)/sqrt(colSums(X^2)))
  rownames(X) = paste("dim", 1:ncol(gmat), sep = "")
  X
}

# Work with rotated matrices and include dimension reduction
computeDelta2 = function(G, means, z0, ndrift=1000){
  require(mvtnorm)
  
  ge = eigen(G)
  G_ge = diag(ge$values[which(ge$values>0)])

  outdat = matrix(NA, nrow=nrow(means), ncol=10)
  
  #Drift expectation
  rb = randomDelta(n = ndrift, gmat = G_ge)
  e_beta = NULL
  e_drift = mean(evolvabilityBeta(G_ge*100, rb)$e)
  
  for(i in 1:nrow(means)){
    z1 = unlist(means[i,]) #Means for population z1
    delta = log(z1)-log(z0) #Difference on log scale from z1 to focal pop (z0)
    rot_delta = t(ge$vectors) %*% delta #Rotate onto the G eigenvectors
    rot_delta = rot_delta[1:ncol(G_ge)] #Reduce dimensions to match reduced G
    scale_delta = rot_delta/sqrt(sum(rot_delta^2)) #Normalize to unit length 
    
    d = rot_delta
    div = mean(abs(d))*100

    e_delta = evolvabilityBeta(G_ge*100, scale_delta)$e
    c_delta = evolvabilityBeta(G_ge*100, scale_delta)$c
    theta = acos(eigen(G_ge)$vectors[,1] %*% scale_delta)*(180/pi)
    
    outdat[i,]=c(div, e_delta, c_delta, theta, evolvabilityMeans(G_ge*100)[c(1:4,7)], e_drift)
  }
  outdat = as.data.frame(outdat)
  colnames(outdat) = c("div", "edelta", "cdelta", "theta", "emean","emin", "emax","cmean","imean","edrift")
  return(outdat)
}

# Add possibility of SE via resampling
computeDelta3 = function(G=out$G, means, z0, SE=FALSE, nSample=10, nFam=out$nFam){
  ge = eigen(G)
  G_ge = diag(ge$values[which(ge$values>0)])
  
  outdat = matrix(NA, nrow=nrow(means), ncol=12)
  for(i in 1:nrow(means)){
    z1 = unlist(means[i,]) #Means for population z1
    delta = log(z1)-log(z0) #Difference on log scale from z1 to focal pop (z0)
    rot_delta = t(ge$vectors) %*% delta #Rotate onto the G eigenvectors
    rot_delta = rot_delta[1:ncol(G_ge)] #Reduce dimensions to match reduced G
    scale_delta = rot_delta/sqrt(sum(rot_delta^2)) #Normalize to unit length 
    
    d = rot_delta
    div = mean(abs(delta))*100
    
    e_delta = evolvabilityBeta(G_ge*100, scale_delta)$e
    c_delta = evolvabilityBeta(G_ge*100, scale_delta)$c
    theta = acos(eigen(G_ge)$vectors[,1] %*% scale_delta)*(180/pi)
    
    outdat[i,1:9]=c(div, e_delta, c_delta, theta, evolvabilityMeans(G_ge*100)[c(1:4,7)])
  }
  
  if(SE){
    require(mvtnorm)
    print("Computing standard errors by resampling")
    se_delta = sc_delta = stheta = matrix(NA, nrow=nSample, ncol=nrow(means))
    
    for(s in 1:nSample){
      sG = round(cov(rmvnorm(nFam, mean=rep(1, ncol(G)), sigma=G)), 10)
      sge = eigen(sG)
      sG_ge = diag(sge$values[which(sge$values>0)])
    
    for(i in 1:nrow(means)){
      z1 = unlist(means[i,]) #Means for population z1
      delta = log(z1)-log(z0) #Difference on log scale from z1 to focal pop (z0)
      rot_delta = t(sge$vectors) %*% delta #Rotate onto the G eigenvectors
      rot_delta = rot_delta[1:ncol(sG_ge)] #Reduce dimensions to match reduced G
      scale_delta = rot_delta/sqrt(sum(rot_delta^2)) #Normalize to unit length 
      
      d = rot_delta
      sdiv = mean(abs(d))*100
      
      se_delta[s,i] = evolvabilityBeta(sG_ge*100, scale_delta)$e/evolvabilityMeans(sG_ge*100)[1]
      sc_delta[s,i] = evolvabilityBeta(sG_ge*100, scale_delta)$c/evolvabilityMeans(sG_ge*100)[4]
      stheta[s,i] = acos(eigen(sG_ge)$vectors[,1] %*% scale_delta)*(180/pi)
    }
    }
      outdat[,10:12]=cbind(apply(se_delta, 2, sd), apply(sc_delta, 2, sd), apply(stheta, 2, sd)) 
    
  }

  
  outdat = as.data.frame(outdat)
  colnames(outdat) = c("div", "edelta", "cdelta", "theta", "emean", "emin", "emax", "cmean", 
                       "imean", "edeltaSE", "cdeltaSE", "thetaSE")
  return(outdat)
}

