estimateDlog = function(means, eV, samples=1000, thin=100){
  
  drop = which(is.na(rowSums(means[-1])))
  if(length(drop)>0){
    means=means[-drop]
    eV=ev[-drop]
  }
  
  #Set prior
  n = ncol(means)
  #alpha.mu <- rep(0, n)
  #alpha.V <- diag(n)*400
  prior<-list(R=list(V=diag(n), nu=n+0.002-1))
  
  mev = eV/(means^2)
  means = apply(means, 2, function(x) log(x)*100)
  mev = melt(mev)$value*10000
  data = as.data.frame(means)
  vars = paste0(colnames(means), collapse=", ")
  vars = paste0("c(",vars,") ~-1+trait")
  
  samples = samples
  thin = thin
  burnin = samples*thin*.5
  nitt = (samples*thin)+burnin
  
  mod<-MCMCglmm(as.formula(noquote(vars)),
                rcov = ~us(trait):units,
                mev = mev,
                data = data, 
                family = rep("gaussian", n), prior = prior, 
                nitt = nitt, burnin = burnin, thin = thin)
  
  modD = matrix(apply(mod$VCV, 2, median)[2:(1+(ncol(means))^2)], nrow=n)/10000
  colnames(modD) = rownames(modD) = colnames(means)
  #modD = meanStdG(modD, colMeans(means))
  return(mod$VCV[, 2:(1+(ncol(means))^2)]/10000)
}
