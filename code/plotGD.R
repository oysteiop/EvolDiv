
######################
Combined

load("data/EVOBASE.RData")
load("data/POPBASE.RData")

#Species present in both databases
gsp=unlist(lapply(EVOBASE, function(x) x$Species))
dsp=unlist(lapply(POPBASE, function(x) x$Species))
both_sp=unique(gsp[which(gsp %in% dsp)])
both_sp

plotGD(species=both_sp[2], gmatrix="mean", dmatrix=1, nbeta=1000,log=T)
plotGD(species=both_sp[2], gmatrix=1, dmatrix=1, log=T)

pdf("figs/GvsDplots.pdf",width=5,height=5)
for(s in 1:length(both_sp)){
  species=both_sp[s]  
  npop=length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
      plotGD(species=species, gmatrix="mean", dmatrix=1,nbeta=1000,log=T)
}    
dev.off()


gmatrix=1
npop=1
species=both_sp[2]
dmatrix=1
nbeta=1000

######################################################################
#### - Plotting function that plots G vs D and returns G and D - #####
######################################################################
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

plotGD=function(species, gmatrix="mean", dmatrix=1, nbeta=1000, betacol="grey",log=T){
  
  npop=length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  
  if(is.numeric(gmatrix) & gmatrix>npop){ 
    stop(paste("Sorry, only", npop, "G matrices available for", species))
         }
  
  if(gmatrix=="mean"){
    gmat=computeMeanG(species)
  }
  
  else{
  gindex=which(unlist(lapply(EVOBASE, function(x) x$Species))==species)[gmatrix]
  gmat=EVOBASE[[gindex]]$G
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
  
  gmatscaled=meanStdG(gmat, means)
  
  dindex=which(unlist(lapply(POPBASE, function(x) x$Species))==species)[dmatrix]
  dmat=POPBASE[[dindex]]$D
  
  m=match(colnames(gmatscaled),colnames(dmat))
  dmat=dmat[m,m]
  
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
    plot(log10(ebeta),log10(dbeta),col=betacol,main=paste(species), las = 1,
         xlab="Evolvability",
         ylab="Population divergence")
         #xlim=c(log10(gmin),log10(gmax)),
         #ylim=c(min(c(log10(dmin),log10(diag(dmat)),log10(dbeta)),na.rm=T),
          #      max(c(log10(dmax),log10(diag(dmat)),log10(dbeta)),na.rm=T)))
    points(log10(diag(gmatscaled)),log10(diag(dmat)),pch=16)
    points(log10(gmax),log10(dmax),col="red", pch=16)
    points(log10(gmin),log10(dmin),col="blue", pch=16)
  }
  
  if(!log){
    plot(ebeta,dbeta,col=betacol,main=paste(species), las=1,
         xlab="Evolvability",
         ylab="Population divergence",
         xlim=c(gmin,gmax),
         ylim=c(min(c(dmin,diag(dmat),dbeta)),max(c(dmax,diag(dmat),dbeta))))
    points(diag(gmatscaled),diag(dmat),pch=16)
    points(gmax,dmax,col="red", pch=16)
    points(gmin,dmin,col="blue", pch=16)
  }
  return(list(D=dmat,G=gmatscaled))
}
