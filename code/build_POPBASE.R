############################################
#### - Building the G matrix database - ####
############################################

library(plyr)
library(reshape2)
library(evolvability)
list.files()
ddat=read.table("data/dmatdata.txt", header=T)
ddat$ID=paste(ddat$reference,ddat$species,ddat$environment, sep="_")

studies=unique(ddat$ID)
studies
length(studies)

Dlist=list()
for(s in 1:length(studies)){
  red=ddat[ddat$ID==studies[[s]],]
  df=dcast(red,population~trait, value.var="mean")
  df=na.omit(df[,-1])
  df=apply(abs(df),2,log) #NB
  
  dmat=cov(df)
  Dlist[[s]]=dmat
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
                      D = signif(Dlist[[i]],4)) 
}

save(POPBASE, file = "data/POPBASE.RData")

