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
for(s in 1:length(studies)){
  
  #Population means
  red=ddat[ddat$ID==studies[[s]],]
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
  
}

names(Dlist)=studies

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
                      popmeans=meanList[[i]],
                      eV = eVlist[[i]],
                      D = signif(Dlist[[i]],4)) 
}

POPBASE[[1]]

save(POPBASE, file = "data/POPBASE.RData")


map1 <- readOGR("C:/data/Political Map.shp")
plot(map1)
points(ddat$lon,ddat$lat,pch=16,col="blue")
