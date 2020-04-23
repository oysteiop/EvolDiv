###########################################################
#### - Computing statistics for scaling of D with G - #####
###########################################################
computeMeanG = function(species){
  sp = unlist(lapply(EVOBASE, function(x) x$Species))
  gMats = lapply(EVOBASE[which(sp==species)], function(x) x$G)
  MeanG = apply(simplify2array(gMats), 1:2, mean)
  return(MeanG)
}

droptraits = function(x){
  drop = which(colSums(x>-Inf, na.rm=T)<max(colSums(x>-Inf, na.rm=T)))
  if(length(drop)>0){
    x = x[-drop,-drop]
  }
  else{
    x = x
  }}

prepareGD = function(species, gmatrix=1, dmatrix=1){
  require(evolvability)
  nG = length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD = length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
  
  if(is.numeric(gmatrix) & gmatrix>nG){ 
    stop(paste("Sorry, only", nG, "G matrices available for", species))
  }
  if(is.numeric(dmatrix) & dmatrix>nD){ 
    stop(paste("Sorry, only", nD, "D matrices available for", species))
  }
  
  gindex = which(unlist(lapply(EVOBASE, function(x) x$Species))==species)[gmatrix]
  gmat = EVOBASE[[gindex]]$G
  gmatname = EVOBASE[[gindex]]$Study_ID
  gmatname = names(EVOBASE)[gindex]
  
  #Remove NAs
  gmat = droptraits(gmat)
  
  #Extract means
  means = EVOBASE[[gindex]]$Means
  means = means[which(names(means) %in% colnames(gmat))]
  
  gmatscaled = meanStdG(gmat, means)
  
  dindex = which(unlist(lapply(POPBASE, function(x) x$Species))==species)[dmatrix]
  dmat = POPBASE[[dindex]]$D
  dmatname = POPBASE[[dindex]]$Study_ID
  dmatname = names(POPBASE)[dindex]

  #Match the matrices
  m = match(colnames(gmatscaled), colnames(dmat))
  dmat = dmat[m,m]
  dmat = droptraits(dmat)
  
  m2 = match(colnames(dmat), colnames(gmatscaled))
  gmatscaled = gmatscaled[m2, m2]
  
  outlist = list(species = species,
                 dmat = dmatname, #Descriptor of D matrix
                 gmat = gmatname, #Descriptor of G matrix
                 D = dmat, #The D matrix
                 G = gmatscaled, #The G matrix
                 nPop = as.numeric(POPBASE[[dindex]]$nPop), #Number of populations included in D
                 nFam = as.numeric(EVOBASE[[gindex]]$nFam)) #Number of families used in estimating G
                 
  return(outlist)
  
}
