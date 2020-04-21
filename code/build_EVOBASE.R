############################################
#### - Building the G matrix database - ####
############################################

rm(list=ls())
library(plyr)
library(evolvability)
source("code/Gvar2G.R")
list.files()

indat = read.table("data/gmatdata.txt", header=T)
sort(unique(indat$studyID))

indat$morph = as.character(indat$morph)
indat$morph[which(is.na(indat$morph))] = "all"
indat$ID = paste(indat$reference, indat$population, indat$environment, indat$morph, sep="_")
studies = unique(indat$ID)
studies

indat$mean[which(indat$e_on_log==1)]=1 #Set mean to 1 for log-scale estimates to get correct mean-scaled G

# Compile G matrices and trait means
Glist = list()
groupList = list()
MeanList = list()
dimList = list()
VpList = list()
isCor = tapply(indat$isCor, indat$ID, mean)>0

for(s in 1:length(studies)){
  red = indat[indat$ID==studies[s],]
  traits = unique(red$traitX)
  #length(traits)
  
  means = NULL
  for(t in 1:length(traits)){
    means[t] = red$mean[which(red$traitX==traits[t] & red$traitY==traits[t])]
  }
  names(means) = traits
  MeanList[[s]] = means
  
  groups = NULL
  for(t in 1:length(traits)){
    groups[t] = as.character(red$traitgroup1[which(red$traitX==traits[t] & red$traitY==traits[t])])
  }
  names(groups) = traits
  groupList[[s]] = groups
  
  dims = NULL
  for(t in 1:length(traits)){
    dims[t] = as.character(red$dimension[which(red$traitX==traits[t] & red$traitY==traits[t])])
  }
  names(dims) = traits
  dimList[[s]] = dims
  
  pvars = NULL
  for(t in 1:length(traits)){
    pvars[t] = red$Vp[which(red$traitX==traits[t] & red$traitY==traits[t])]
  }
  names(pvars) = traits
  VpList[[s]] = pvars
  
  G = matrix(NA, nrow = length(traits), ncol = length(traits))
  for(i in 1:length(traits)){
    for(j in 1:length(traits)){
      w = max(which(red$traitX==traits[i] & red$traitY==traits[j]),
              which(red$traitX==traits[j] & red$traitY==traits[i]))
      G[i,j] = red$Va[w]
    }
  }
  colnames(G) = rownames(G) = traits
  Glist[[s]] = G
  
  if(isCor[studies[[s]]]){  
    G = matrix(NA, nrow = length(traits), ncol = length(traits))
    for(i in 1:length(traits)){
      for(j in 1:length(traits)){
        w = max(which(red$traitX==traits[i] & red$traitY==traits[j]),
                which(red$traitX==traits[j] & red$traitY==traits[i]))
        G[i,j] = ifelse(max(c(red$corA[w],red$h2[w]),na.rm=T)>-2,max(c(red$corA[w],red$h2[w]),na.rm=T),NA)
      }
    }
    
    G = Gvar2G(G, VpList[[s]])
    colnames(G) = rownames(G) = traits
    Glist[[s]] = G
    
  }
}

# Preparing metadata
names(indat)
metadata = ddply(indat, .(ID), summarize,
                 Family = family[1],
                 Species = species[1],
                 Population = population[1],
                 Lat = mean(lat, na.rm=T),
                 Lon = mean(lon, na.rm=T),
                 Environment = environment[1],
                 nfam = median(nfam, na.rm=T))
head(metadata)

#### Building the database list ####
EVOBASE = list()
for(i in 1:length(studies)){
  EVOBASE[[i]] = list(Study_ID = paste(metadata$ID[metadata$ID==studies[i]]), 
                      Family = paste(metadata$Family[metadata$ID==studies[i]]),
                      Species = paste(metadata$Species[metadata$ID==studies[i]]), 
                      Population = paste(metadata$Population[metadata$ID==studies[i]]),
                      LatLon = c(Lat = metadata$Lat[metadata$ID==studies[i]],
                                 Lon = metadata$Lon[metadata$ID==studies[i]]),
                      Environment = paste(metadata$Environment[metadata$ID==studies[i]]),
                      G = signif(Glist[[i]],4), 
                      Groups=groupList[[i]],
                      Dims=dimList[[i]],
                      Means = MeanList[[i]],
                      Vp = VpList[[i]])
}

sp_names = unlist(lapply(EVOBASE, function(x) x$Species))
pop_names = unlist(lapply(EVOBASE, function(x) x$Population))
titles = paste0(sp_names,": ", pop_names)

for(i in 1:length(titles)){
  if(duplicated(titles)[i]){
    titles[i] = paste(titles[i], "II", sep=" ")
  }
  if(duplicated(titles)[i]){
    titles[i] = paste0(titles[i], "I")
  }
  if(duplicated(titles)[i]){
    titles[i] = paste0(titles[i], "I")
  }
}

titles

names(EVOBASE) = gsub("_", " ", titles)

EVOBASE = EVOBASE[order(titles)]

save(EVOBASE, file = "data/EVOBASE.RData")

View(EVOBASE)
