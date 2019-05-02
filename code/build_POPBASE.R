############################################
#### - Building the D matrix database - ####
############################################

rm(list=ls())

library(plyr)
library(reshape2)
library(evolvability)
library(rgdal)

list.files()
ddat=read.table("data/dmatdata.txt", header=T)
ddat$ID=paste(ddat$reference,ddat$species,ddat$environment, sep="_")

pops_per_study=tapply(ddat$population, ddat$ID, function(x) length(unique(x)))

ddat=ddat[-which(ddat$ID%in%names(which(pops_per_study<2))),]
studies=unique(ddat$ID)
studies
length(studies)

meanList=list()
Dlist=list()
eVlist=list()
distmatList=list()
for(s in 1:length(studies)){
  
  #Population means
  red=ddat[ddat$ID==studies[[s]],]
  red$population=factor(red$population)
  df=dcast(red,population~trait, value.var="mean")
  meanList[[s]]=df
  
  #D matrix
  df=na.omit(df[,-1])
  df=apply(abs(df),2,log) #NB
  dmat=cov(df)
  Dlist[[s]]=dmat
  
  #Error variances
  df2=dcast(red,population~trait, value.var="vp")
  #df2=na.omit(df2[,-1])
  n=dcast(red,population~trait, value.var="n")
  df2[,-1]=df2[,-1]/n[,-1]
  eVlist[[s]]=df2
  
  #Distance matrices
  require(geosphere)
  distdat=cbind(tapply(red$lon, red$population, mean, na.rm=T),
            tapply(red$lat, red$population, mean, na.rm=T))
  distmat=matrix(NA,nrow(distdat),nrow(distdat))
  for(i in 1:nrow(distdat)){
    distmat[,i] = round((distVincentyEllipsoid(c(distdat[i,1],distdat[i,2]),cbind(distdat[,1],distdat[,2]))),2)
  }
  rownames(distmat)=colnames(distmat)=rownames(distdat)
  
  distmatList[[s]]=distmat
  
}

names(Dlist)=studies
meanList[3]
eVlist[[3]]

######Estimate error-corrected D matrices ####

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

####
means=meanList[[18]]
eV=eVlist[[18]]

drop=which(is.na(rowSums(means[-1])))
if(length(drop)>0){
means=means[-drop]
eV=ev[-drop]
}

#Set prior
n = ncol(means)-1
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1))

means[,2:ncol(means)]=apply(means[,2:ncol(means)], 2, function(x) x*100)

mev=melt(eV[,-1])$value*10000

data=means[,-1]

vars=paste0(colnames(means)[-1], collapse=", ")
vars=paste0("c(",vars,") ~-1+trait")

mod<-MCMCglmm(as.formula(noquote(vars)),
              #random = ~us(trait):OutcropSystem,
              rcov = ~us(trait):units,
              mev = mev,
              data = means, 
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

#summary(mod)
plot(mod$VCV[,2])

modD=matrix(apply(mod$VCV, 2, median)[2:(1+(ncol(means)-1)^2)]/10000,nrow=n)
colnames(modD) = rownames(modD) = colnames(means)[-1]
round(modD,3)

#cov(means[,-1]/100)

#plot(c(modD),c(cov(means[,-1]/100)))
#lines(0:100,0:100)
  
#adjDlist[s]=modD


#Preparing metadata
metadata = ddply(ddat, .(ID), summarize,
                 Family = family[1],
                 Species = species[1],
                 nPop = length(unique(population)))

metadata

#### - Building the database list - ####
POPBASE = list()
for(i in 1:length(studies)){
  POPBASE[[i]] = list(Study_ID = paste(metadata$ID[metadata$ID==studies[i]]), 
                      Family = paste(metadata$Family[metadata$ID==studies[i]]),
                      Species = paste(metadata$Species[metadata$ID==studies[i]]), 
                      nPop = paste(metadata$nPop[metadata$ID==studies[i]]),
                      distmat=round(distmatList[[i]],0),
                      popmeans=meanList[[i]],
                      eV = eVlist[[i]],
                      D = signif(Dlist[[i]],4)) 
}

POPBASE[[1]]

save(POPBASE, file = "data/POPBASE.RData")


map1 <- readOGR("C:/data/Political Map.shp")
plot(map1)
points(ddat$lon,ddat$lat,pch=16,col="blue")
