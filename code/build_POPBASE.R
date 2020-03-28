############################################
#### - Building the D matrix database - ####
############################################

rm(list=ls())

library(plyr)
library(reshape2)
library(evolvability)
library(rgdal)

list.files()
ddat = read.table("data/dmatdata.txt", header=T)
ddat$ID = paste(ddat$reference, ddat$species, ddat$environment, sep="_")
ddat$upop = paste(ddat$species, ddat$population, sep="_")
ddat$utrait = paste(ddat$species, ddat$trait, sep="_")

# Summary stats
sum(ddat$mean>0)
length(unique(ddat$species))
length(unique(ddat$upop))
length(unique(ddat$utrait))
tapply(ddat$upop, ddat$environment, function(x) length(unique(x)))

names(ddat)

pops_per_study = tapply(ddat$population, ddat$ID, function(x) length(unique(x)))

ddat = ddat[-which(ddat$ID%in%names(which(pops_per_study<2))),]
studies = unique(ddat$ID)
studies
length(studies)

meanList=list()
Dlist=list()
eVlist=list()
distmatList=list()

for(s in 1:length(studies)){
  
  #Population means
  red = ddat[ddat$ID==studies[[s]],]
  red$population = factor(red$population)
  df = dcast(red, population~trait, value.var="mean")
  meanList[[s]] = df
  
  #D matrix
  df = na.omit(df[,-1])
  df = apply(abs(df),2,log) #NB
  dmat = cov(df)
  Dlist[[s]] = dmat
  
  #Error variances
  df2 = dcast(red, population~trait, value.var="vp")
  #df2=na.omit(df2[,-1])
  n = dcast(red, population~trait, value.var="n")
  df2[,-1] = df2[,-1]/n[,-1]
  eVlist[[s]] = df2
  
  #Distance matrices
  require(geosphere)
  distdat = cbind(tapply(red$lon, red$population, mean, na.rm=T),
                  tapply(red$lat, red$population, mean, na.rm=T))
  distmat = matrix(NA, nrow(distdat), nrow(distdat))
  for(i in 1:nrow(distdat)){
    distmat[,i] = round((distVincentyEllipsoid(c(distdat[i,1], distdat[i,2]), cbind(distdat[,1], distdat[,2]))), 2)
  }
  rownames(distmat) = colnames(distmat) = rownames(distdat)
  
  distmatList[[s]] = distmat
}

names(Dlist) = studies
meanList[39]
eVlist[[3]]

# Preparing metadata
metadata = ddply(ddat, .(ID), summarize,
                 Family=family[1],
                 Species=species[1],
                 nPop=length(unique(population)))
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

length(POPBASE)
POPBASE[[39]]

save(POPBASE, file = "data/POPBASE.RData")

load(file = "data/POPBASE.RData")

map1 <- readOGR("C:/data/Political Map.shp")
plot(map1)
points(ddat$lon,ddat$lat,pch=16,col="blue")

maxdists = cbind(print(unlist(lapply(POPBASE, function(x)x$Study_ID))),
      print(unlist(lapply(POPBASE, function(x) max(x$distmat)/1000))))
View(maxdists)

maxdists[which(maxdists[,1]=="Barrett_&_Shore_1987_Turnera_ulmifolia_greenhouse"),2]=188 #Ca. distance from map in paper
maxdists[which(maxdists[,1]=="Podolsky_et_al._1997_Clarkia_dudleyana_greenhouse"),2]=11.3 #Ca. distance from map in paper
maxdists[which(maxdists[,1]=="Carter_&_Murdy_1986_Talinum_mengesii_greenhouse"),2]=214 #Ca. distance from map in paper

maxdists[which(maxdists[,1]=="Carr_&_Fenster_1994_Mimulus_guttatus_greenhouse"),2]=10 #Description in paper
maxdists[which(maxdists[,1]=="Fenster_&_Carr_1997_Mimulus_guttatus_greenhouse"),2]=10 #Description in Carr & Fenster 1994
maxdists[which(maxdists[,1]=="Caruso_2000_Ipomopsis_aggregata_field"),2]=3.5 #Ca. distance author description/map in paper
maxdists[which(maxdists[,1]=="Caruso_2001_Ipomopsis_aggregata_field"),2]=3 #Ca. distance ("a few km") from author description
  
maxdists[which(maxdists[,1]=="Andersson_Crepis_Crepis_tectorum_greenhouse"),2]=2500 #Ca. distance author description, without the Canadian pops
maxdists[which(maxdists[,1]=="Billington_et_al_1988_Holcus_lanatus_greenhouse"),2]=.5 #Ca. distance from description in paper ("adjacent fields")
maxdists[which(maxdists[,1]=="Mcgoey_&_Stinchcombe_2018_Ambrosia_artemisiifolia_common_garden"),2]=540 #Ca. distance for North American pops from map in paper. Around 180 km for French pops
maxdists[which(maxdists[,1]=="Colautti_&_Barrett_2011_Lythrum_salicaria_greenhouse"),2]=1200 #Ca. distance from coordinate range in paper

maxdists[which(maxdists[,1]=="Campbell_et_al_2018_Ipomopsis_tenuituba_field"),2]=0.0943 #Author description in Excel file
maxdists[which(maxdists[,1]=="Campbell_et_al_2018_Ipomopsis_aggregata_x_tenuituba_field"),2]=0.96 #Author description in Excel file
maxdists[which(maxdists[,1]=="Campbell_et_al_2018_Ipomopsis_aggregata_field"),2]=0.199 #Author description in Excel file

write.csv(maxdists, file="data/maxdists.csv")
