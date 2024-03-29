############################################
#### - Building the D matrix database - ####
############################################

rm(list=ls())

library(plyr)
library(reshape2)
library(evolvability)
library(rgdal)

list.files()
#ddat = read.csv2("data/dmatdata.csv", dec=".")
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

meanList = list()
Dlist = list()
eVlist = list()
distmatList = list()

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
  if(sum(is.na(eVlist[[s]]))>0){
    df2 = dcast(red, population~trait, value.var="se")
    eVlist[[s]] = df2[,-1]^2
  }
  
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
meanList[41]
eVlist[[41]]

# Preparing metadata
metadata = ddply(ddat, .(ID), summarize,
                 Family=family[1],
                 Species=species[1],
                 nPop=length(unique(population)),
                 Env = unique(environment))
metadata

#### - Building the database list - ####
POPBASE = list()
for(i in 1:length(studies)){
  POPBASE[[i]] = list(Study_ID = paste(metadata$ID[metadata$ID==studies[i]]), 
                      Family = paste(metadata$Family[metadata$ID==studies[i]]),
                      Species = paste(metadata$Species[metadata$ID==studies[i]]), 
                      nPop = paste(metadata$nPop[metadata$ID==studies[i]]),
                      distmat=round(distmatList[[i]],0),
                      Env = paste(metadata$Env[metadata$ID==studies[i]]),
                      popmeans=meanList[[i]],
                      eV = eVlist[[i]],
                      D = signif(Dlist[[i]],4)) 
}

length(POPBASE)

sp_names = unlist(lapply(POPBASE, function(x) x$Species))
study_names = unlist(lapply(POPBASE, function(x) x$Study_ID))

for(i in 1:length(study_names)){
  study_names[i]=gsub(pattern = sp_names[i], rep="", x=study_names[i])
  study_names[i]=gsub(pattern = "__", rep="_", x=study_names[i])
}
study_names
titles = paste0(sp_names,": ", study_names)
titles

names(POPBASE) = gsub("_", " ", titles)

POPBASE = POPBASE[order(titles)]

save(POPBASE, file = "data/POPBASE.RData")

View(POPBASE)

load(file = "data/POPBASE.RData")

maxdists = cbind(print(unlist(lapply(POPBASE, function(x)x$Study_ID))),
      print(unlist(lapply(POPBASE, function(x) max(x$distmat)/1000))))


maxdists[which(maxdists[,1]=="Barrett_and_Shore_1987_Turnera_ulmifolia_greenhouse"),2]=188 #Ca. distance from map in paper
maxdists[which(maxdists[,1]=="Podolsky_et_al._1997_Clarkia_dudleyana_greenhouse"),2]=11.3 #Ca. distance from map in paper
maxdists[which(maxdists[,1]=="Carter_and_Murdy_1986_Talinum_mengesii_greenhouse"),2]=214 #Ca. distance from map in paper
maxdists[which(maxdists[,1]=="Carter_and_Murdy_1986_Talinum_teretifolium_greenhouse"),2]=110 #Ca. distance from map in paper

maxdists[which(maxdists[,1]=="Carr_and_Fenster_1994_Mimulus_guttatus_greenhouse"),2]=10 #Description in paper
maxdists[which(maxdists[,1]=="Fenster_and_Carr_1997_Mimulus_guttatus_greenhouse"),2]=10 #Description in Carr & Fenster 1994
maxdists[which(maxdists[,1]=="Caruso_2000_Ipomopsis_aggregata_field"),2]=3.5 #Ca. distance author description/map in paper
maxdists[which(maxdists[,1]=="Caruso_2001_Ipomopsis_aggregata_field"),2]=3 #Ca. distance ("a few km") from author description
  
maxdists[which(maxdists[,1]=="Andersson_Crepis_Crepis_tectorum_greenhouse"),2]=2500 #Ca. distance author description, without the Canadian pops
maxdists[which(maxdists[,1]=="Billington_et_al_1988_Holcus_lanatus_greenhouse"),2]=.5 #Ca. distance from description in paper ("adjacent fields")
maxdists[which(maxdists[,1]=="Mcgoey_and_Stinchcombe_2018_Ambrosia_artemisiifolia_common_garden"),2]=540 #Ca. distance for North American pops from map in paper. Around 180 km for French pops
maxdists[which(maxdists[,1]=="Colautti_and_Barrett_2011_Lythrum_salicaria_greenhouse"),2]=1200 #Ca. distance from coordinate range in paper

maxdists[which(maxdists[,1]=="Campbell_et_al_2018_Ipomopsis_tenuituba_field"),2]=0.0943 #Author description in Excel file
maxdists[which(maxdists[,1]=="Campbell_et_al_2018_Ipomopsis_aggregata_x_tenuituba_field"),2]=0.96 #Author description in Excel file
maxdists[which(maxdists[,1]=="Campbell_et_al_2018_Ipomopsis_aggregata_field"),2]=0.199 #Author description in Excel file

maxdists=as.data.frame(maxdists)
maxdists$ID = gsub("_", " ", titles)

View(maxdists)

write.csv(maxdists, file="data/maxdists.csv")

