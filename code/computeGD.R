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
  
  outlist=list(dmat=dmatname,
               gmat=gmatname,
               D=dmat,
               G=gmatscaled,
               nPop=as.numeric(POPBASE[[dindex]]$nPop),
               betaG=beta_g,
               r2G=r2_g,
               betaD=beta_d,
               r2D=r2_d,
               iG=evolvabilityMeans(gmatscaled)[7],
               nBetaG=length(mg$residuals),
               nBetaD=length(md$residuals))
               
  #Plot
  #xmin=log10(min(c(var_g_g, var_g_d), na.rm=T))
  #xmax=log10(max(c(var_g_g, var_g_d), na.rm=T))
  #ymin=log10(min(c(var_d_g, var_d_d), na.rm=T))
  #ymax=log10(max(c(var_d_g, var_d_d), na.rm=T))
  #plot(log10(var_g_g), log10(var_d_g), xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  #points(log10(var_g_d), log10(var_d_d), pch=16)
  #legend("bottomright", c("G eigenvectors", "D eigenvectors"), pch=c(1,16))
  
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

computeGD(species = both_sp[1], gmatrix = 1, dmatrix = 1)

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
    for(d in 1:nD[s])
      res=computeGD(species = both_sp[s], gmatrix = g, dmatrix = d)[c(1:2, 5:12)]
      reslist[length(reslist)+1]=as.data.frame(unlist(res)[c(1:10)])
  }
}

reslist[[1]]

resmat=as.data.frame(matrix(NA, ncol=10, nrow=29))
for(i in 1:29){
  resmat[i,1:2]=as.character(reslist[[i]][1:2])
  resmat[i,3:8]=round(as.numeric(as.character(reslist[[i]][3:8])),2)
  resmat[i,9:10]=as.numeric(as.character(reslist[[i]][9:10]))
 }

colnames(resmat)=c("D", "G", "npop", "betaG", "r2G", "betaD", "r2D","i_mean","nBetaG","nBetaD")
resmat=resmat[c(2,1,3:10)]
head(resmat)
View(resmat)


hist(resmat$betaG)
hist(resmat$r2G[resmat$nBetaG>2])
