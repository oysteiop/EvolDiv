species=both_sp[5]
gmatrix=1
dmatrix=1

###########################################################
#### - Computing statistics for scaling of D with G - #####
###########################################################
computeMeanG=function(species){
  sp=unlist(lapply(EVOBASE, function(x) x$Species))
  gMats=lapply(EVOBASE[which(sp==species)], function(x) x$G)
  MeanG=apply(simplify2array(gMats), 1:2, mean)
  return(MeanG)
}

droptraits=function(x){
  drop=which(colSums(x>-Inf,na.rm=T)<max(colSums(x>-Inf,na.rm=T)))
  if(length(drop)>0){
    x=x[-drop,-drop]
  }
  else{
    x=x
  }}

computeGD=function(species, gmatrix=1, dmatrix=1){
  require(evolvability)
  nG=length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD=length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
  
  if(is.numeric(gmatrix) & gmatrix>nG){ 
    stop(paste("Sorry, only", nG, "G matrices available for", species))
  }
  if(is.numeric(dmatrix) & dmatrix>nD){ 
    stop(paste("Sorry, only", nD, "D matrices available for", species))
  }
  
  gindex=which(unlist(lapply(EVOBASE, function(x) x$Species))==species)[gmatrix]
  gmat=EVOBASE[[gindex]]$G
  gmatname=EVOBASE[[gindex]]$Study_ID
  
  #Remove NAs
  gmat=droptraits(gmat)
  
  #Extract means
  means=EVOBASE[[gindex]]$Means
  means=means[which(names(means) %in% colnames(gmat))]
  
  gmatscaled=meanStdG(gmat, means)
  
  dindex=which(unlist(lapply(POPBASE, function(x) x$Species))==species)[dmatrix]
  dmat=POPBASE[[dindex]]$D
  dmatname=POPBASE[[dindex]]$Study_ID
  
  #Match the matrices
  m=match(colnames(gmatscaled),colnames(dmat))
  dmat=dmat[m,m]
  dmat=droptraits(dmat)
  
  m2=match(colnames(dmat), colnames(gmatscaled))
  gmatscaled=gmatscaled[m2,m2]
  
  #Compute eigenvectors etc.
  g_ev=eigen(gmatscaled)$vectors
  var_g_g=evolvabilityBeta(gmatscaled, Beta = g_ev)$e
  var_d_g=evolvabilityBeta(dmat, Beta = g_ev)$e
  
  d_ev=eigen(dmat)$vectors
  var_g_d=evolvabilityBeta(gmatscaled, Beta = d_ev)$e
  var_d_d=evolvabilityBeta(dmat, Beta = d_ev)$e

  #Compute summary stats
  mg=lm(log(var_d_g)~log(var_g_g))
  beta_g=summary(mg)$coef[2,1]
  r2_g=summary(mg)$r.squared
  
  md=lm(log(var_d_d)~log(var_g_d))
  beta_d=summary(md)$coef[2,1]
  r2_d=summary(md)$r.squared
  
  mt=lm(log(diag(dmat))~log(diag(gmatscaled)))
  beta_t=summary(mt)$coef[2,1]
  r2_t=summary(mt)$r.squared
  
  outlist=list(dmat = dmatname, #Descriptor of D matrix
               gmat = gmatname, #Descriptor of G matrix
               D = dmat, #The D matrix
               G = gmatscaled, #The G matrix
               nPop = as.numeric(POPBASE[[dindex]]$nPop), #Number of populations included in D
               theta=acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi), #Angle between leading eigenvectors of G and D
               betaG = beta_g, #Slope of D on G for eigenvectors of G
               r2G = r2_g, #r^2 for betaG
               betaD = beta_d, #Slope for D on G for eigenvectors of D
               r2D = r2_d, #r^2 for betaD
               betaT = beta_t, #Slope for traits
               r2T = r2_t, #r^2 for traits
               iG = evolvabilityMeans(gmatscaled)[7], #Hansen-Houle integration index for G
               emean = evolvabilityMeans(gmatscaled)[1]*100, #Mean evolvability
               dmean = evolvabilityMeans(dmat)[1]*100, #Mean divergence
               nBetaG = length(mg$residuals), #Number of positive eigenvalues compared for G eigenvectors
               nBetaD = length(md$residuals), #Number of positive eigenvalues compared for D eigenvectors
               nBetaT = length(mt$residuals)) #Number of traits compared
               
    return(outlist)
  
  }

######################################
######################################
load("data/EVOBASE.RData")
load("data/POPBASE.RData")

#Species present in both databases
gsp=unlist(lapply(EVOBASE, function(x) x$Species))
dsp=unlist(lapply(POPBASE, function(x) x$Species))
both_sp=unique(gsp[which(gsp %in% dsp)])
both_sp

POPBASE=POPBASE[-c(14:15)] #Dropping second M. guttatus study

computeGD(species = both_sp[2], gmatrix = 1, dmatrix = 1)

nG=NULL
nD=NULL
for(s in c(1:length(both_sp))){
  species=both_sp[s]  
  nG[s]=length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD[s]=length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
}
cbind(both_sp, nG, nD)

reslist=list()
for(s in 1:length(both_sp)){
  for(g in 1:nG[s]){
    for(d in 1:nD[s]){
      res=computeGD(species = both_sp[s], gmatrix = g, dmatrix = d)[c(1:2, 5:18)]
      reslist[length(reslist)+1]=as.data.frame(unlist(res)[c(1:16)])
    }
  }
}
length(reslist)
reslist[[20]]

resmat=as.data.frame(matrix(NA, ncol=16, nrow=length(reslist)))
for(i in 1:length(reslist)){
  resmat[i,1:2]=as.character(reslist[[i]][1:2])
  resmat[i,3:13]=round(as.numeric(as.character(reslist[[i]][3:13])),2)
  resmat[i,14:16]=as.numeric(as.character(reslist[[i]][14:16]))
 }

colnames(resmat)=c("D", "G", "npop", "theta", "betaG", "r2G", "betaD", "r2D","betaT","r2T","i_mean","e_mean","d_mean","nBetaG","nBetaD","nBetaT")

for(i in 1:nrow(resmat)){
  if(resmat$theta[i]>=90){
    resmat$theta[i] = 180 - resmat$theta[i]
  }
}

names(resmat)
resmat=resmat[,c(2,1,16,3, 12,13,5,6)]
head(resmat)
View(resmat)

resmat=resmat[resmat$npop>2,]
resmat=resmat[resmat$nBetaG>2,]
resmat=resmat[resmat$nBetaD>2,]
resmat=resmat[resmat$nBetaT>2,]

mean(resmat$betaG)
hist(resmat$r2G)

plot(resmat$betaT, resmat$betaG)
plot(resmat$r2T, resmat$r2G)

plot(resmat$theta, resmat$betaG)

plot(resmat$npop, resmat$betaG)
plot(resmat$npop, resmat$r2G, cex=resmat$nBetaG*.5)

plot(resmat$i_mean, resmat$r2G, cex=resmat$nBetaG*.5, xlim=c(0,1))
