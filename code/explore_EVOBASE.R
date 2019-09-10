##########################################
#### - Extracting data from EVOBASE - ####
##########################################
load("data/EVOBASE.RData")

#Extracting a single species
print(unlist(lapply(EVOBASE, function(x)x$Study_ID)))
EVOBASE[[34]]

lapply(EVOBASE, function(x) x$Dims)

#Extracting just the G matrices
GG=lapply(EVOBASE,function(x) x$G)

#Function to remove NAs from G matrices
droptraits=function(x){
  drop=which(colSums(x>-Inf,na.rm=T)<max(colSums(x>-Inf,na.rm=T)))
  if(length(drop)>0){
    x=x[-drop,-drop]
  }
  else{
    x=x
  }}

GG=lapply(GG,droptraits)

#Extracting the means
mm=lapply(EVOBASE,function(x) x$Means)

#Remove trait means for traits not in G
for(i in 1:length(mm)){
  mm[[i]]=mm[[i]][which(names(mm[[i]]) %in% colnames(GG[[i]]))]
}

#Mean-scale all the G matrices
MG=list()
for(i in 1:length(GG)){
  MG[[i]]=GG[[i]]/tcrossprod(mm[[i]],mm[[i]])*100
}

#Running some function on each mean-scaled G matrix
elist=lapply(MG,evolvabilityMeans)
names(elist)=lapply(EVOBASE, function(x) x$Study_ID)
elist

