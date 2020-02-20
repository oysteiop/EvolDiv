gmatrix=1
npop=1
species=both_sp[12]
dmatrix=1
nbeta=1000

######################################################################
#### - Plotting function that plots G vs D and returns G and D - #####
######################################################################
computeMeanG=function(species){
  sp=unlist(lapply(EVOBASE, function(x) x$Species))
  gMats=lapply(EVOBASE[which(sp==species)], function(x) x$G)
  MeanG=apply(simplify2array(gMats), 1:2, mean, na.rm=T) #NB, with na.rm=T will use whatever data is available
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

plotGD=function(species, gmatrix="mean", dmatrix=1, nbeta=1000, betacol="grey", log=T){
  require(evolvability)
  nG=length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD=length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
  
  if(is.numeric(gmatrix) & gmatrix>nG){ 
    stop(paste("Sorry, only", nG, "G matrices available for", species))
  }
  if(is.numeric(dmatrix) & dmatrix>nD){ 
    stop(paste("Sorry, only", nD, "D matrices available for", species))
  }
  
  if(gmatrix=="mean"){
    gmat=computeMeanG(species)
    gmatname="Mean"
  }
  
  else{
    gindex=which(unlist(lapply(EVOBASE, function(x) x$Species))==species)[gmatrix]
    gmat=EVOBASE[[gindex]]$G
    gmatname=EVOBASE[[gindex]]$Study_ID
  }
  
  #Remove NAs
  gmat=droptraits(gmat)
  
  #Extract means
  if(gmatrix=="mean"){
    means=lapply(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)],function(x) x$Means)
    means=apply(simplify2array(means),1,mean)
  }
  
  else{
    means=EVOBASE[[gindex]]$Means
  }
  means=means[which(names(means) %in% colnames(gmat))]
  
  gmatscaled=meanStdG(gmat, means)*100
  
  dindex=which(unlist(lapply(POPBASE, function(x) x$Species))==species)[dmatrix]
  dmat=POPBASE[[dindex]]$D*100
  dmatname=POPBASE[[dindex]]$Study_ID
  
  #Match the matrices
  m=match(colnames(gmatscaled),colnames(dmat))
  dmat=dmat[m,m]
  dmat=droptraits(dmat)
  
  m2=match(colnames(dmat), colnames(gmatscaled))
  gmatscaled=gmatscaled[m2,m2]
  
  
  #Compute eigenvectors etc.
  first_ev=eigen(gmatscaled)$vectors[,1]
  gmax=evolvabilityBeta(gmatscaled, Beta = first_ev)$e
  dmax=evolvabilityBeta(dmat, Beta = first_ev)$e
  
  last_ev=eigen(gmatscaled)$vectors[,nrow(gmatscaled)]
  gmin=evolvabilityBeta(gmatscaled, Beta = last_ev)$e
  dmin=evolvabilityBeta(dmat, Beta = last_ev)$e
  
  betas=randomBeta(nbeta,nrow(gmatscaled))
  ebeta=evolvabilityBeta(gmatscaled, betas)$e
  dbeta=evolvabilityBeta(dmat, betas)$e
  
  if(log){
    xminvals=na.omit(c(log10(gmin), log10(diag(gmatscaled)), log10(ebeta)))
    xmaxvals=na.omit(c(log10(gmax), log10(diag(gmatscaled)), log10(ebeta)))
    yminvals=na.omit(c(log10(dmin), log10(diag(dmat)), log10(dbeta)))
    ymaxvals=na.omit(c(log10(dmax), log10(diag(dmat)), log10(dbeta)))
    
    
        plot(log10(ebeta), log10(dbeta), col=betacol, main=paste(species), las = 1,
         xlab="Evolvability (log %)",
         ylab="Population divergence (log %)",
          xlim=c(min(xminvals[xminvals>-Inf]),max(xmaxvals[xmaxvals>-Inf])),
          ylim=c(min(yminvals[yminvals>-Inf]),max(ymaxvals[ymaxvals>-Inf])))
    
    points(log10(diag(gmatscaled)),log10(diag(dmat)),pch=16)
    points(log10(gmax),log10(dmax),col="red", pch=16)
    points(log10(gmin),log10(dmin),col="blue", pch=16)
  }
  
  if(!log){
    plot(ebeta, dbeta, col=betacol, main=paste(species), las=1,
         xlab="Evolvability (%)",
         ylab="Population divergence (%)",
         xlim=c(min(c(gmin,diag(gmatscaled),ebeta)),max(c(gmax,diag(gmatscaled),ebeta))),
         ylim=c(min(c(dmin,diag(dmat),dbeta)),max(c(dmax,diag(dmat),dbeta))))
    points(diag(gmatscaled),diag(dmat),pch=16)
    points(gmax,dmax,col="red", pch=16)
    points(gmin,dmin,col="blue", pch=16)
  }
  
  legend("topleft", c(paste("G:",gmatname), paste("D:",dmatname)), bty="n",cex=.7)
  legend("bottomright",c("gmax","gmin"),pch=16,col=c("red","blue"))
  return(list(D=dmat,G=gmatscaled))
}



######################################
#### - Code for plotting G vs D - ####
######################################
load("data/EVOBASE.RData")
load("data/POPBASE.RData")

#Species present in both databases
gsp=unlist(lapply(EVOBASE, function(x) x$Species))
dsp=unlist(lapply(POPBASE, function(x) x$Species))
both_sp=unique(gsp[which(gsp %in% dsp)])
both_sp

x11()
par(mfrow=c(1,2))
plotGD(species=both_sp[11], gmatrix=1, dmatrix=2, nbeta=1000, log=F, betacol="white")
plotGD(species=both_sp[11], gmatrix=1, dmatrix=3, nbeta=1000, log=F, betacol="white")

#plotGD(species=both_sp[11], gmatrix="mean", dmatrix=4, nbeta=1000,log=F)


POPBASE = POPBASE[-c(14:15)] #Dropping second M. guttatus study

pdf("figs/GvsDplots.pdf", width=8, height=4.5)
for(s in c(1:length(both_sp))){
  species = both_sp[s]  
  nG = length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD = length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
    for(g in 1:nG){  
      for(d in 1:nD){  
        par(mfrow=c(1,2))
        plotGD(species=species, gmatrix=g, dmatrix=d, nbeta=1000, log=F)
        plotGD(species=species, gmatrix=g, dmatrix=d, nbeta=1000, log=T)
      }
    }
  }
dev.off()


