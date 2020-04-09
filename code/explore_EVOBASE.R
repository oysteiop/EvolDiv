##########################################
#### - Extracting data from EVOBASE - ####
##########################################

rm(list=ls())
library(evolvability)

# Function to remove NAs from G matrices
droptraits=function(x){
  drop=which(colSums(x>-Inf,na.rm=T)<max(colSums(x>-Inf,na.rm=T)))
  if(length(drop)>0){
    x=x[-drop,-drop]
  }
  else{
    x=x
  }}

load("data/EVOBASE.RData")

# List studies in the database
print(unlist(lapply(EVOBASE, function(x) x$Study_ID)))

# Extracting a single species
s=52
View(EVOBASE[[s]])

EVOBASE[[s]]$G
EVOBASE[[s]]$Means

# Mean-scale the G matrix
mg = meanStdG(EVOBASE[[s]]$G, EVOBASE[[s]]$Means)*100
mg

evolvabilityMeans(mg)
eigen(mg)
signif(cov2cor(mg), 3)

# Extracting just all the G matrices
GG = lapply(EVOBASE, function(x) x$G)


GG = lapply(GG, droptraits)

# Extracting the means
mm = lapply(EVOBASE, function(x) x$Means)

# Remove trait means for traits not in G
for(i in 1:length(mm)){
  mm[[i]]=mm[[i]][which(names(mm[[i]]) %in% colnames(GG[[i]]))]
}

# Mean-scale all the G matrices
MG=list()
for(i in 1:length(GG)){
  MG[[i]]=GG[[i]]/tcrossprod(mm[[i]],mm[[i]])*100
}

# Running some function on each mean-scaled G matrix
elist = lapply(MG, evolvabilityMeans)
names(elist) = lapply(EVOBASE, function(x) x$Study_ID)
elist

#### Reducing to linear floral traits ####
GG = lapply(EVOBASE, function(x){sel = which(x$Groups=="floral" & x$Dims=="linear") 
                                   x$G[sel, sel]})

sel = which(unlist(lapply(GG, function(x) dim(as.matrix(x))[1]))>1)
GG = GG[sel]
GG = lapply(GG, droptraits)

# Extracting the means
mm = lapply(EVOBASE, function(x){sel = which(x$Groups=="floral" & x$Dims=="linear") 
x$Means[sel]})

mm = mm[sel]

# Remove trait means for traits not in G
for(i in 1:length(mm)){
  mm[[i]]=mm[[i]][which(names(mm[[i]]) %in% colnames(GG[[i]]))]
}

# Mean-scale all the G matrices
MG=list()
for(i in 1:length(GG)){
  MG[[i]]=GG[[i]]/tcrossprod(mm[[i]],mm[[i]])*100
}

names(MG) = lapply(EVOBASE, function(x) x$Study_ID)[sel]
MG

evolvabilityBeta(MG[37][[1]], Beta=c(1,0,0))
cov2cor(MG[34][[1]])

elist = lapply(MG, evolvabilityMeans)
elist

