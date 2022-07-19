############################################################
#### Analyses of evolvability along divergence vectors #####
############################################################

rm(list=ls())

add2deltaList=function(pop=POPBASE[[s]]$popmeans[,1]){
  data.frame(species = out$species, g = out$gmat, ntraits = ncol(out$G),
             d = out$dmat, pop = pop, 
             dims = paste0(substr(sort(unique(EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))])), 1, 4), collapse="+"),
             ndims = length(unique(EVOBASE[[as.character(out$g)]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))])),
             traitgroups = paste0(substr(sort(unique(EVOBASE[[out$g]]$Groups[match(colnames(out$G), names(EVOBASE[[out$g]]$Groups))])), 1, 3), collapse="+"),
             emean = outdat$emean,
             emin = outdat$emin,
             emax = outdat$emax,
             cmean = outdat$cmean,
             div = outdat$div, edelta = outdat$edelta, cdelta = outdat$cdelta,
             theta = outdat$theta, edrift = outdat$edrift, row.names=NULL)
}

library(plyr)
library(reshape2)
library(MCMCglmm)
library(evolvability)
source("code/prepareGD.R")
source("code/computeDelta.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")

# Species present in both databases
gsp = unlist(lapply(EVOBASE, function(x) x$Species))
dsp = unlist(lapply(POPBASE, function(x) x$Species))
names(dsp)
both_sp = unique(gsp[which(gsp %in% dsp)])
both_sp

deltaList = list()

#### Lobelia ####

# CERA G matrix, D = Caruso 2003 field
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means
colnames(means)==names(z0)

#outdat = computeDelta(G=out$G, means=means, z0=z0)
#evolvabilityMeans(out$G)
outdat = computeDelta2(G=out$G, means=means, z0=z0)
#outdat2 = computeDelta3(G=out$G, means=means, z0=z0, SE=T, nSample=100)

deltaList[[length(deltaList)+1]] = add2deltaList()

# CERA G matrix, D=Caruso 2012 
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Means
s=out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1,3:6)] #Caruso 2012
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

# Krumm G matrix, D=Caruso 2003
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 3)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means #Caruso 2003
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

# Krumm G matrix, D=Caruso 2012 
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1,3:6)] #Caruso 2012
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

#### Aquilegia ####

# G = QFP1, D = Herlihy and Eckert 2007
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(means[2,])
means = means[-2,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-2,1])

# G = QFP1, D = Bartkowska et al. 2018
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# G = QLL3, D = Herlihy and Eckert 2007
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(means[3,]) #Second G mat, first D mat
means = means[-3,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-3,1])

# G = QLL3, D = Bartkowska et al. 2018
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means #Second D mat
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

#### Ipomopsis ####
# D = Caruso 2000
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0=EVOBASE[[out$gmat]]$Means[c(2:5)]
colnames(means)==names(z0)

outdat=computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# D = Caruso 2001
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 4)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(2,3,5)]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# D = Campbell 2018
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1:3)]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# D = Campbell 2019 unpubl.
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1:3, 5)]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

#### Brassica ####
gg=1
for(gg in 1:5){
  out = prepareGD(species = "Brassica_cretica", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s = out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(means[gg,])
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(G=out$G, means=means, z0=z0)
  #outdat2 = computeDelta3(G=out$G, means=means, z0=z0, SE=T, nSample=100)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Turnera ####
out = prepareGD(species = "Turnera_ulmifolia", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

#### Clarkia ####
out = prepareGD(species = "Clarkia_dudleyana", gmatrix = 2, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(c(means[4,])) #Second G-matrix
means = means[-4,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-4,1])

#### Spergularia ####
for(gg in 1:4){
  out = prepareGD(species = "Spergularia_marina", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s=out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Solanum ####
for(gg in 1:3){
  out = prepareGD(species = "Solanum_carolinense", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s=out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  
  colnames(means)==names(z0)
  
  outdat=computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Fragaria ####

# G = Heterophrodite
out = prepareGD(species = "Fragaria_virginiana", gmatrix = 3, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = apply(means, 2, mean)
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# G = Female
out = prepareGD(species = "Fragaria_virginiana", gmatrix = 4, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = apply(means, 2, mean)
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

##### Mimulus guttatus ####
gmats = c(5,7)
gg=1
for(gg in 1:2){
out = prepareGD(species = "Mimulus_guttatus", gmatrix = gmats[gg], dmatrix = 1)
out$gmat
out$dmat

# Means
names(POPBASE)
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(c(means[gg,])) #First G-matrix
means = means[-gg,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)
outdat2 = computeDelta3(out$G, means, z0, SE=T, nSample=100)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}


##### Mimulus micranthus ####
for(gg in 1:2){
  out = prepareGD(species = "Mimulus_micranthus", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Meangg
  s = out$dmat
  means = abs(POPBASE[[s]]$popmeans[,-1])
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

##### Nigella degenii ####
gg=1
for(gg in 1:2){
  out = prepareGD(species = "Nigella_degenii", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s = out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Lobelia ####
for(gg in 1:2){
  out = prepareGD(species = "Lobelia_siphilitica", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s = out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}


#### Plotting deltaList ####

#save(deltaList, file="deltaList.RData")

load(file="deltaList.RData")

# Add data from orginal analyses
load(file="analyses/andersson_crepis/deltaDF.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/dalechampia/deltaDF_Tulum_CR.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/dalechampia/deltaDF_Tulum_MX.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/dalechampia/deltaDF_Tovar_MX.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/walter2018/deltaDF_Dune.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/walter2018/deltaDF_Head.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/walter2018/deltaDF_Table.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/walter2018/deltaDF_Wood.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/puentes2016/deltaDF_SPIT.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/puentes2016/deltaDF_STUC.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/puentes2016/deltaDF_STUS.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/puentes2016/deltaDF_VIS.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/colautti/deltaDF.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/delph_Silene/deltaDF.RData")
deltaList[[length(deltaList)+1]]=deltaDF

save(deltaList, file="deltaListFull.RData")

# Compile dataframe
deltaDat = rbind.fill(deltaList)
dim(deltaDat)
tail(deltaDat)

# Keep angles between 0 and 90
for(i in 1: nrow(deltaDat)){
  if(deltaDat$theta[i]>90){
    deltaDat$theta[i]=180-deltaDat$theta[i]
  }
}

save(deltaDat, file="deltaDat.RData")

#### END OF FINAL ANALYSES ####
