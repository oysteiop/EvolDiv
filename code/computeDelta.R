computeDelta = function(G, means, z0){
outdat = matrix(NA, nrow=nrow(means), ncol=4)
for(i in 1:nrow(means)){
  z1 = unlist(means[i,])
  delta = log(z1)-log(z0)
  scale_delta = delta/sqrt(sum(delta^2)) 
  d = delta
  div = mean(abs(d))*100
  
  e_delta = evolvabilityBeta(G*100, scale_delta)$e
  c_delta = evolvabilityBeta(G*100, scale_delta)$c
  theta = acos(t(eigen(G)$vectors[,1]) %*% scale_delta)*(180/pi)
  
  outdat[i,]=c(div, e_delta, c_delta, theta)
  }
return(outdat)
}
