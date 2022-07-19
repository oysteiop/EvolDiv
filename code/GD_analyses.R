#######################################
#### Analyses of G and D matrices #####
#######################################

rm(list=ls())

add2gdList=function(){
  data.frame(species = out$species, g = out$gmat, ntraits = ncol(out$G),
             nPop = out$nPop, nFam = out$nFam,
             dims = paste0(substr(sort(unique(EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))])), 1, 4), collapse="+"),
             ndims = length(unique(EVOBASE[[as.character(out$g)]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))])),
             traitgroups = paste0(substr(sort(unique(EVOBASE[[out$g]]$Groups[match(colnames(out$G), names(EVOBASE[[out$g]]$Groups))])), 1, 3), collapse="+"),
             emean = vals$evolvabilityMeans[1],
             emin = vals$evolvabilityMeans[2],
             emax = vals$evolvabilityMeans[3],
             cmean = vals$evolvabilityMeans[4],
             imean = vals$evolvabilityMeans[7],
             d = out$dmat, 
             dmean = vals$divergenceMeans[1],
             betaT = vals$res[1,3], betaT_SE = vals$res[1,5], r2T = vals$res[1,6],
             betaT_cond = vals$res[2,3], betaT_cond_SE = vals$res[2,5], r2T_cond = vals$res[2,6],
             betaG = vals$res[3,3], betaG_SE = vals$res[3,5], r2G = vals$res[3,6],
             betaD = vals$res[4,3], betaD_SE = vals$res[4,5], r2D = vals$res[4,6],
             betaD_cond = vals$res[5,3], betaD_cond_SE = vals$res[5,5], r2D_cond = vals$res[5,6],
             betaP = vals$res[6,3], betaP_SE = vals$res[6,5], r2P = vals$res[6,6],
             betaP_cond = vals$res[7,3], betaP_cond_SE = vals$res[7,5], r2P_cond = vals$res[7,6],
             r2All = vals$res[8,6],
             theta = vals$theta, row.names = NULL)
}

library(reshape2)
library(plyr)
library(MCMCglmm)
library(lme4)
library(evolvability)
source("code/prepareGD.R")
source("code/computeGD.R")
#source("code/computeGDorg.R")
#source("code/estimateD.R")
source("code/estimateDlog.R")
source("code/alignMat.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")
load("data/PMATBASE.RData")

# Species present in both databases
gsp = unlist(lapply(EVOBASE, function(x) x$Species))
dsp = unlist(lapply(POPBASE, function(x) x$Species))
both_sp = unique(gsp[which(gsp %in% dsp)])
both_sp

estD = FALSE #Re-estimating the error-corrected D matrices?
thin = 100 # Thinning interval for the models estimating error-corrected D

fixD = FALSE #Hold D fixed when assessing uncertainty
nSample = 10 #Number of resamples for SE

linearonly = F

gdList = list()

# Start of Analyses ####

#### Lobelia ####

# G = CERA, D = Caruso 2012
out = prepareGD(species="Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)$G
glist[[2]] = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

# P matrix
colnames(out$G)
names(PMATBASE)
plist = list()
ma = match(colnames(out$G), names(PMATBASE[["Lobelia siphilitica: CERA"]]$Means))

plist[[1]] = meanStdG(PMATBASE[["Lobelia siphilitica: CERA"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: CERA"]]$Means[ma])
plist[[2]] = meanStdG(PMATBASE[["Lobelia siphilitica: Krumm"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Krumm"]]$Means[ma])
plist[[3]] = meanStdG(PMATBASE[["Lobelia siphilitica: Reichelt"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Reichelt"]]$Means[ma])
MeanP = apply(simplify2array(plist), 1:2, mean)

#EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))]

vals = computeGD(out$G, modDpost*100, MeanP, species="Lobelia siphilitica", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]]=add2gdList()

# G = CERA, D = Caruso 2003
out = prepareGD(species="Lobelia_siphilitica", gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)$G
glist[[2]] = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)  #Fails with longer thin but looks OK with 100
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Lobelia_Caruso2003.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2003.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

# P matrix
colnames(out$G)
names(PMATBASE)
plist = list()
ma = match(colnames(out$G), names(PMATBASE[["Lobelia siphilitica: CERA"]]$Means))
plist[[1]] = meanStdG(PMATBASE[["Lobelia siphilitica: CERA"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: CERA"]]$Means[ma])
plist[[2]] = meanStdG(PMATBASE[["Lobelia siphilitica: Krumm"]]$P[ma], 
                      PMATBASE[["Lobelia siphilitica: Krumm"]]$Means[ma])
plist[[3]] = meanStdG(PMATBASE[["Lobelia siphilitica: Reichelt"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Reichelt"]]$Means[ma])
MeanP = apply(simplify2array(plist), 1:2, mean)

#vals = computeGD(out$G, out$D, MeanP, species="Lobelia siphilitica", SE=T, fixD=T, nSample=nSample, plot=F)
vals = computeGD(out$G, modDpost, MeanP, species="Lobelia siphilitica", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]]=add2gdList()

# G = Krumm, D = Caruso 2012
out = prepareGD(species="Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

# P matrix
colnames(out$G)
names(PMATBASE)
plist = list()
ma = match(colnames(out$G), names(PMATBASE[["Lobelia siphilitica: CERA"]]$Means))
plist[[1]] = meanStdG(PMATBASE[["Lobelia siphilitica: CERA"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: CERA"]]$Means[ma])
plist[[2]] = meanStdG(PMATBASE[["Lobelia siphilitica: Krumm"]]$P[ma], 
                      PMATBASE[["Lobelia siphilitica: Krumm"]]$Means[ma])
plist[[3]] = meanStdG(PMATBASE[["Lobelia siphilitica: Reichelt"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Reichelt"]]$Means[ma])
MeanP = apply(simplify2array(plist), 1:2, mean)

vals = computeGD(out$G, modDpost, MeanP, species="Lobelia siphilitica", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

# G = Krumm, D = Caruso 2003
out = prepareGD(species="Lobelia_siphilitica", gmatrix = 2, dmatrix = 3)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2003.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

# P matrix
colnames(out$G)
names(PMATBASE)
plist = list()
ma = match(colnames(out$G), names(PMATBASE[["Lobelia siphilitica: CERA"]]$Means))
plist[[1]] = meanStdG(PMATBASE[["Lobelia siphilitica: CERA"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: CERA"]]$Means[ma])
plist[[2]] = meanStdG(PMATBASE[["Lobelia siphilitica: Krumm"]]$P[ma], 
                      PMATBASE[["Lobelia siphilitica: Krumm"]]$Means[ma])
plist[[3]] = meanStdG(PMATBASE[["Lobelia siphilitica: Reichelt"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Reichelt"]]$Means[ma])
MeanP = apply(simplify2array(plist), 1:2, mean)

vals = computeGD(out$G, modDpost, MeanP, species="Lobelia siphilitica", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#### Aquilegia ####

# G = QFP1, D = Bartkowska 2018
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)$G
glist[[2]] = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Aquilegia_Bartkowska2018.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_Bartkowska2018.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

# P matrix
indat = read.csv("data/eckert/Hierstruct_floral_morphology.csv")
pops = unique(indat$OutcropSystem)

colnames(out$G)

plist=list()
for(i in 1:length(pops)){
  pop = pops[i]
  df = na.omit(subset(indat[indat$OutcropSystem==pop,], 
            select = c("PistilLength.mm", "SpurLength.mm", "Herkogamy.mm")))
  plist[[i]] = meanStdG(cov(df), colMeans(df))
}
MeanP = apply(simplify2array(plist),1:2,mean)
colnames(MeanP) = rownames(MeanP) = colnames(out$G)

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G),2)
#signif(cov2cor(out$D),2)

vals = computeGD(out$G, modDpost, MeanP, species="Aquilegia canadensis", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#G = QLL3, D = Bartkowska 2018
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_Bartkowska2018.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G),2)
#signif(cov2cor(out$D),2)

vals = computeGD(out$G, modDpost, MeanP, species="Aquilegia canadensis", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#G = QFP1, D = Herlihy and Eckert 2007
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)$G
glist[[2]] = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Aquilegia_HerlihyEckert2007.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_HerlihyEckert2007.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G),2)
#signif(cov2cor(out$D),2)

vals = computeGD(out$G, modDpost, MeanP, species="Aquilegia canadensis", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#G = QLL3, D = Herlihy and Eckert 2007
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_HerlihyEckert2007.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G),2)
#signif(cov2cor(out$D),2)

vals = computeGD(out$G, modDpost, MeanP, species="Aquilegia canadensis", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#### Brassica ####
out = prepareGD(species="Brassica_cretica", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Mean G matrix
glist = list()
for(i in 1:5){
  glist[[i]] = prepareGD(species="Brassica_cretica", gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Brassica.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Brassica.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

# Loop
SEvals = c(T,T,T,F,T)
gg=5
for(gg in c(1:5)){
out = prepareGD(species="Brassica_cretica", gmatrix = gg, dmatrix = 1)

load(file="analyses/adj_Dmats/Brassica.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100 #Single G

#vals = computeGD(out$G, out$D, species="Brassica cretica", SE=F, fixD=T,  nSample=nSample, plot="e", xmin=-3, ymin=-2)
vals = computeGD(out$G, modDpost, species="Brassica cretica", linearonly=linearonly, SE=SEvals[gg], fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()
}

#### Spergularia ####
out = prepareGD(species="Spergularia_marina", gmatrix = 4, dmatrix = 1)
out$G

# Mean G matrix
glist = list()
for(i in 1:4){
  glist[[i]] = prepareGD(species = "Spergularia_marina", gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD*100)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Spergularia.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Spergularia.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

# P matrix
names(PMATBASE)

plist = list()
ma = match(colnames(out$G), names(PMATBASE[["Spergularia marina: ASC"]]$Means))

plist[[1]] = meanStdG(PMATBASE[["Spergularia marina: ASC"]]$P[ma, ma], 
                      PMATBASE[["Spergularia marina: ASC"]]$Means[ma])
plist[[2]] = meanStdG(PMATBASE[["Spergularia marina: COP"]]$P[ma, ma], 
                      PMATBASE[["Spergularia marina: COP"]]$Means[ma])
plist[[3]] = meanStdG(PMATBASE[["Spergularia marina: MSH"]]$P[ma, ma], 
                      PMATBASE[["Spergularia marina: MSH"]]$Means[ma])
plist[[4]] = meanStdG(PMATBASE[["Spergularia marina: SMB"]]$P[ma, ma], 
                      PMATBASE[["Spergularia marina: SMB"]]$Means[ma])
MeanP = apply(simplify2array(plist), 1:2, mean)
colnames(out$G)

# Loop
gg=1
for(gg in 1:4){
out = prepareGD(species="Spergularia_marina", gmatrix = gg, dmatrix = 1)

load(file="analyses/adj_Dmats/Spergularia.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100 #Single G

vals = computeGD(out$G, modDpost, MeanP, species="Spergularia marina", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, ymin=-4, plot=F)

gdList[[length(gdList)+1]] = add2gdList()
}

#### Solanum ####
out = prepareGD(species = "Solanum_carolinense", gmatrix = 1, dmatrix = 1)
out$gmat

# Mean G matrix
glist = list()
for(i in 1:3){
  glist[[i]] = prepareGD(species = "Solanum_carolinense", gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin) #Fails with thin 1000 but looks OK with 100
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Solanum.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Solanum.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

# Loop
for(gg in 1:3){
out = prepareGD(species = "Solanum_carolinense", gmatrix = gg, dmatrix = 1)

load(file="analyses/adj_Dmats/Solanum.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100 #Single G

EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))]

vals = computeGD(out$G, modDpost, species="Solanum carolinense", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()
}

#### Clarkia ####
out = prepareGD(species = "Clarkia_dudleyana", gmatrix = 2, dmatrix = 1)
out$gmat

# Mean G matrix
glist = list()
for(i in 1:2){
  glist[[i]] = prepareGD(species = "Clarkia_dudleyana", gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 4)*100
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Clarkia.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Clarkia.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)
EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))]

vals = computeGD(out$G, modDpost, species="Clarkia dudleyana", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#### Ipomopsis ####

# D = Campbell 2019 unpubl.
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 2)
out$dmat

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Ipomopsis_Campbell2019.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Campbell2019.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100
#out$D = out$D*100

# P matrix
colnames(out$G)
names(PMATBASE)
ma = match(colnames(out$G), names(PMATBASE[["Ipomopsis aggregata: Vera Falls"]]$Means))
MeanP = meanStdG(PMATBASE[["Ipomopsis aggregata: Vera Falls"]]$P[ma, ma], PMATBASE[["Ipomopsis aggregata: Vera Falls"]]$Means[ma])

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)

EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))]
#out$G = out$G[-4, -4]
#out$D = out$D[-4, -4]
#MeanP = MeanP[-4, -4]

vals = computeGD(out$G, modDpost, MeanP, species="Ipomopsis aggregata", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

# D = Caruso 2000
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 3)
out$dmat

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Ipomopsis_Caruso2000.RData")
}

# P matrix
colnames(out$G)
names(PMATBASE)
ma = match(colnames(out$G), names(PMATBASE[[1]]$Means))

plist = list()
for(i in 1:7){
  plist[[i]] = meanStdG(PMATBASE[[i]]$P[ma, ma], PMATBASE[[i]]$Means[ma])
}
MeanP = apply(simplify2array(plist), 1:2, mean)
  
# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2000.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)

vals = computeGD(out$G, modDpost, MeanP, species="Ipomopsis aggregata", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

# D = Caruso 2001
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 4)
out$dmat

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Ipomopsis_Caruso2001.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2001.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100
#out$D = out$D*100

# P matrix
colnames(out$G)
names(PMATBASE)
ma = match(colnames(out$G), names(PMATBASE[[1]]$Means))

plist = list()
for(i in 1:7){
  plist[[i]] = meanStdG(PMATBASE[[i]]$P[ma, ma], PMATBASE[[i]]$Means[ma])
}
MeanP = apply(simplify2array(plist), 1:2, mean)

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)

vals = computeGD(out$G, modDpost, MeanP, species="Ipomopsis aggregata", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#### Turnera ####
out = prepareGD(species = "Turnera_ulmifolia", gmatrix = 1, dmatrix = 1)
out$gmat

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Turnera.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Turnera.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)
EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))]

vals = computeGD(out$G, modDpost, species="Turnera ulmifolia", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#### Fragaria ####

#NB: G-matrices are for females and hermaphrodites separately, D is for overall means

# Hermaphrodite G
out = prepareGD(species = "Fragaria_virginiana", gmatrix = 3, dmatrix = 1)
out$gmat

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Fragaria.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Fragaria.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)
EVOBASE[[out$g]]$Dims[match(colnames(out$G), names(EVOBASE[[out$g]]$Dims))]

vals = computeGD(out$G, modDpost, species="Fragaria virginiana", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

# Female G
out = prepareGD(species = "Fragaria_virginiana", gmatrix = 4, dmatrix = 1)
out$gmat

# Estimate error-corrected D matrix
if(estD){
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modDpost = estimateDlog(means, eV, thin=thin)
modD = matrix(apply(modDpost, 2, median), nrow=length(means))
colnames(modD) = rownames(modD) = colnames(means)

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

save(modDpost, file="analyses/adj_Dmats/Fragaria_Female.RData")
}

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Fragaria_Female.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = out$G*100

#evolvabilityMeans(out$G)
#evolvabilityMeans(out$D)
#signif(cov2cor(out$G), 2)
#signif(cov2cor(out$D), 2)

vals = computeGD(out$G, modDpost, species="Fragaria virginiana", linearonly=linearonly, SE=T, fixD=fixD, nSample=nSample, plot=F)

gdList[[length(gdList)+1]] = add2gdList()

#### dgList output ####

save(gdList, file="gdListSE.RData")

load(file="gdListSE.RData")

# Add data from orginal analyses
load(file="analyses/dalechampia/gdDF_Tulum_MX.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/dalechampia/gdDF_Tulum_CR.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/dalechampia/gdDF_Tovar_MX.RData")
gdList[[length(gdList)+1]] = gdDF

load(file="analyses/andersson_crepis/gdDF.RData")
gdList[[length(gdList)+1]] = gdDF

load(file="analyses/walter2018/gdDF_Dune.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/walter2018/gdDF_Head.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/walter2018/gdDF_Table.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/walter2018/gdDF_Wood.RData")
gdList[[length(gdList)+1]] = gdDF

load(file="analyses/puentes2016/gdDF_SPIT.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/puentes2016/gdDF_STUC.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/puentes2016/gdDF_STUS.RData")
gdList[[length(gdList)+1]] = gdDF
load(file="analyses/puentes2016/gdDF_VIS.RData")
gdList[[length(gdList)+1]] = gdDF

load(file="analyses/colautti/gdDF.RData")
gdList[[length(gdList)+1]] = gdDF

gdDat = rbind.fill(gdList)
dim(gdDat)
head(gdDat)

for(i in 1: nrow(gdDat)){
  if(gdDat$theta[i]>90){
    gdDat$theta[i]=180-gdDat$theta[i]
  }
}

save(gdDat, file="gdDatSE.RData")

#### END OF FINAL ANALYSES ####
