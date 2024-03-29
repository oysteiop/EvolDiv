##########################################
#### - Extracting data from EVOBASE - ####
##########################################

rm(list=ls())
library(evolvability)
library(maps)

# Function to remove NAs from G matrices
droptraits = function(x){
  drop = which(colSums(x>-Inf, na.rm=T)<max(colSums(x>-Inf, na.rm=T)))
  if(length(drop)>0){
    x = x[-drop,-drop]
  }
  else{
    x = x
  }}

load("data/EVOBASE.RData")

# Map of study populations
coords = na.omit(matrix(unlist(lapply(EVOBASE, function(x) x$LatLon)), ncol=2, byrow=T))

plot=F
if(plot){
x11(height=5, width=7)
map("world", fill=T, col=rgb(0,.8,0,.2), mar=c(2,4,1,1))
map.axes(bty="l", las=1)
points(coords[,2], coords[,1], pch=16, col="darkblue")
}

# List studies in the database
names(EVOBASE)

unique(unlist(lapply(EVOBASE, function(x) x$Species)))

# Matrix of trait groups
ss = lapply(EVOBASE, function(x) x$Subgroups)
ns = lapply(ss, length)
studyVec = rep(names(ss), times=ns)

df = data.frame(study=studyVec, subgroup=as.factor(unlist(ss)))
df = droplevels(df)

sdf = dcast(df, study~subgroup, fun.aggregate = length)

# Extracting a single species
s = "Raphanus raphanistrum: Binghamton III"
View(EVOBASE[[s]])

names(EVOBASE[[s]])
signif(EVOBASE[[s]]$G, 2)
EVOBASE[[s]]$Means
EVOBASE[[s]]$Subgroups

signif(cov2cor(EVOBASE[[s]]$G), 2)

# Mean-scale the G matrix
colnames(EVOBASE[[s]]$G)==names(EVOBASE[[s]]$Means)
mg = meanStdG(EVOBASE[[s]]$G, EVOBASE[[s]]$Means)*100
mg = droptraits(mg)
mg

evolvabilityMeans(mg)
eigen(mg)$values
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

drop = which(lapply(GG, function(x) sum(is.na(x)))>0)
GG = GG[-drop]

# Extracting the means
mm = lapply(EVOBASE, function(x){sel = which(x$Groups=="floral" & x$Dims=="linear") 
                                 x$Means[sel]})
mm = mm[sel]
mm = mm[-drop]

names(GG) == names(mm)

# Remove trait means for traits not in G
for(i in 1:length(mm)){
  mm[[i]]=mm[[i]][which(names(mm[[i]]) %in% colnames(GG[[i]]))]
}

# Mean-scale all the G matrices
MG=list()
for(i in 1:length(GG)){
  MG[[i]]=GG[[i]]/tcrossprod(mm[[i]], mm[[i]])*100
}

names(MG) = lapply(EVOBASE, function(x) x$Study_ID)[sel][-drop]
MG

