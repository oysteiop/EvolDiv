############################################
#### - Building the P matrix database - ####
############################################

library(plyr)
library(evolvability)
list.files()

indat = read.table("data/pmatdata.txt", header=T)
indat$morph = as.character(indat$morph)
indat$morph[which(is.na(indat$morph))] = "all"
indat$ID = paste(indat$reference,indat$population,indat$environment,indat$morph, sep="_")
studies = unique(indat$ID)
studies

#Compile P matrices and trait means
Plist = list()
MeanList = list()
VpList = list()
isCor = tapply(indat$isCor, indat$ID, mean)>0

for(s in 1:length(studies)){
  red = indat[indat$ID==studies[s],]
  traits = unique(red$traitX)
  length(traits)
  
  means = NULL
  for(t in 1:length(traits)){
    means[t] = red$mean[which(red$traitX==traits[t] & red$traitY==traits[t])]
  }
  names(means) = traits
  MeanList[[s]] = means
  
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
  Plist[[s]] = G
  
  
  if(isCor[studies[[s]]]){  
    G = matrix(NA, nrow = length(traits), ncol = length(traits))
    for(i in 1:length(traits)){
      for(j in 1:length(traits)){
        w = max(which(red$traitX==traits[i] & red$traitY==traits[j]),
                which(red$traitX==traits[j] & red$traitY==traits[i]))
        G[i,j] = ifelse(max(c(red$corA[w],red$h2[w]),na.rm=T)>-2,max(c(red$corA[w],red$h2[w]),na.rm=T),NA)
      }
    }
    
    Gvar2G=function(Gvar,Vp){
      vars=diag(Gvar)*Vp
      G=Gvar*sqrt(tcrossprod(vars))
      diag(G)=vars
      return(G)
    }
    
    G=Gvar2G(G,VpList[[s]])
    colnames(G) = rownames(G) = traits
    Plist[[s]] = G
    
  }
}


#Preparing metadata
metadata = ddply(indat, .(ID), summarize,
                 Family = family[1],
                 Species = species[1],
                 Population = population[1],
                 Lat = mean(lat, na.rm=T),
                 Lon = mean(lon, na.rm=T))
metadata

#### - Building the database list - ####
PMATBASE = list()
for(i in 1:length(studies)){
  PMATBASE[[i]] = list(Study_ID = paste(metadata$ID[metadata$ID==studies[i]]), 
                      Family = paste(metadata$Family[metadata$ID==studies[i]]),
                      Species = paste(metadata$Species[metadata$ID==studies[i]]), 
                      Population = paste(metadata$Population[metadata$ID==studies[i]]),
                      LatLon = c(Lat = metadata$Lat[metadata$ID==studies[i]],
                                 Lon = metadata$Lon[metadata$ID==studies[i]]),
                      P = signif(Plist[[i]],4), 
                      Means = MeanList[[i]],
                      Vp = VpList[[i]])
}

sp_names = unlist(lapply(PMATBASE, function(x) x$Species))
pop_names = unlist(lapply(PMATBASE, function(x) x$Population))
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

names(PMATBASE) = gsub("_", " ", titles)

PMATBASE = PMATBASE[order(titles)]

save(PMATBASE, file = "data/PMATBASE.RData")

View(PMATBASE)
