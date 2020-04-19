#######################################
#### Analyses of G and D matrices #####
#######################################

rm(list=ls())
library(reshape2)
library(plyr)
library(MCMCglmm)
library(evolvability)
source("code/computeGD.R")
source("code/plot_GD.R")
source("code/estimateD.R")
source("code/estimateDlog.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")

# Species present in both databases
gsp = unlist(lapply(EVOBASE, function(x) x$Species))
dsp = unlist(lapply(POPBASE, function(x) x$Species))
both_sp = unique(gsp[which(gsp %in% dsp)])
both_sp

#POPBASE = POPBASE[-c(14:15)] #Dropping second M. guttatus study

gdList = list()

#### Lobelia ####

# G = CERA, D = Caruso 2012
out = computeGD(species="Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = computeGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)$G
glist[[2]] = computeGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

#modD = estimateD(means, eV, thin=500)
modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

vals = plot_GD(out$G, out$D, species="Lobelia siphilitica", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Lobelia_siphilitica", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)


# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

# Plot
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Lobelia siphilitica: Caruso 2003 D", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

cvals=NULL
for(i in 1:ncol(out$G)){
  b=rep(0,ncol(out$G))
  b[i]=1
  cvals[i] = evolvabilityBeta(out$G, b)$c
}
points(log10(cvals), log10(diag(out$D)), pch=16, col="red")

# G = CERA, D = Caruso 2003
out = computeGD(species="Lobelia_siphilitica", gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = computeGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)$G
glist[[2]] = computeGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=100)

save(modD, file="analyses/adj_Dmats/Lobelia_Caruso2003.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2003.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

vals = plot_GD(out$G, out$D, species="Lobelia siphilitica", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Lobelia_siphilitica", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)
# G = Krumm, D = Caruso 2012
out = computeGD(species="Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

vals = plot_GD(out$G, out$D, species="Lobelia siphilitica", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Lobelia_siphilitica", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

# G = Krumm, D = Caruso 2003
out = computeGD(species="Lobelia_siphilitica", gmatrix = 2, dmatrix = 3)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2003.RData")

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Individual G

vals = plot_GD(out$G, out$D, species="Lobelia siphilitica", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Lobelia_siphilitica", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#### Aquilegia ####

#G = QFP1, D = Bartkowska 2018
out = computeGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = computeGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)$G
glist[[2]] = computeGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=500)
save(modD, file="analyses/adj_Dmats/Aquilegia_Bartkowska2018.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_Bartkowska2018.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G),2)
signif(cov2cor(out$D),2)

vals = plot_GD(out$G, out$D, species="Aquilegia canadensis")

gdList[[length(gdList)+1]] = data.frame(sp = "Aquilegia_canadensis", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#G = QLL3, D = Bartkowska 2018
out = computeGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_Bartkowska2018.RData")

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G),2)
signif(cov2cor(out$D),2)

vals = plot_GD(out$G, out$D, species="Aquilegia canadensis")

gdList[[length(gdList)+1]] = data.frame(sp = "Aquilegia_canadensis", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)
#G = QFP1, D = Herlihy and Eckert 2007
out = computeGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = computeGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)$G
glist[[2]] = computeGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=100)
save(modD, file="analyses/adj_Dmats/Aquilegia_HerlihyEckert2007.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_HerlihyEckert2007.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G),2)
signif(cov2cor(out$D),2)

vals = plot_GD(out$G, out$D, species="Aquilegia canadensis", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Aquilegia_canadensis", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#G = QLL3, D = Herlihy and Eckert 2007
out = computeGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Aquilegia_HerlihyEckert2007.RData")

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G),2)
signif(cov2cor(out$D),2)

vals = plot_GD(out$G, out$D, species="Aquilegia canadensis")

gdList[[length(gdList)+1]] = data.frame(sp = "Aquilegia_canadensis", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#### Brassica ####

out = computeGD(species="Brassica_cretica", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Mean G matrix
glist = list()
for(i in 1:5){
  glist[[i]] = computeGD(species="Brassica_cretica", gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Brassica.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Brassica.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

# Loop
for(gg in 1:5){
out = computeGD(species="Brassica_cretica", gmatrix = gg, dmatrix = 1)

load(file="analyses/adj_Dmats/Brassica.RData")
out$D = modD*100
out$G = out$G*100 #Single D

vals = plot_GD(out$G, out$D, species="Brassica cretica", xmin=-1, ymin=-2, plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Brassica_cretica", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)
}

#### Spergularia ####
out = computeGD(species="Spergularia_marina", gmatrix = 1, dmatrix = 1)
out$gmat

# Mean G matrix
glist = list()
for(i in 1:4){
  glist[[i]] = computeGD(species = both_sp[14], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

#modD = estimateD(means, eV, thin=100)
modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Spergularia.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Spergularia.RData")

round(modD, 3)
eigen(modD*100)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

# Loop
for(gg in 1:4){
out = computeGD(species="Spergularia_marina", gmatrix = gg, dmatrix = 1)

load(file="analyses/adj_Dmats/Spergularia.RData")
out$D = modD*100
out$G = out$G*100 #Single G

vals = plot_GD(out$G, out$D, species="Spergularia marina", ymin=-4, plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Spergularia_marina", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)
}

#### Solanum ####
out = computeGD(species = "Solanum_carolinense", gmatrix = 1, dmatrix = 1)
out$gmat

# Mean G matrix
glist = list()
for(i in 1:3){
  glist[[i]] = computeGD(species = both_sp[13], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s=out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=100)

save(modD, file="analyses/adj_Dmats/Solanum.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Solanum.RData")

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

# Loop
for(gg in 1:3){
  
out = computeGD(species = "Solanum_carolinense", gmatrix = gg, dmatrix = 1)

load(file="analyses/adj_Dmats/Solanum.RData")
out$D = modD*100
out$G = out$G*100 #Single G

vals = plot_GD(out$G, out$D, species="Solanum carolinense", xmin=-1, ymin=-2)

gdList[[length(gdList)+1]] = data.frame(sp = "Solanum_carolinense", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)
}

#### Clarkia ####
out = computeGD(species = "Clarkia_dudleyana", gmatrix = 2, dmatrix = 1)
out$gmat

# Mean G matrix
glist = list()
for(i in 1:2){
  glist[[i]] = computeGD(species = both_sp[3], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Clarkia.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Clarkia.RData")

round(modD, 4)*100
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100

#out$G = gmat*100 #Mean G
out$G = out$G*100 #Single G

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

vals = plot_GD(out$G, out$D, species="Clarkia dudleyana", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Clarkia_dudleyana", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#### Ipomopsis ####

# D = Campbell 2019 unpubl.
out = computeGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 2)
out$dmat

s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Ipomopsis_Campbell2019.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Campbell2019.RData")

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = out$G*100
#out$D = out$D*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

vals = plot_GD(out$G, out$D, species="Ipomopsis aggregata", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Ipomopsis_aggregata", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)


# D = Caruso 2000
out = computeGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 3)
out$dmat

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Ipomopsis_Caruso2000.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2000.RData")

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = out$G*100
#out$D = out$D*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

vals = plot_GD(out$G, out$D, species="Ipomopsis aggregata", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Ipomopsis_aggregata", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

# D = Caruso 2001
out = computeGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 4)
out$dmat

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Ipomopsis_Caruso2001.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2001.RData")

round(modD, 3)
eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = out$G*100
#out$D = out$D*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

vals = plot_GD(out$G, out$D, species="Ipomopsis aggregata", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Ipomopsis_aggregata", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#### Turnera ####
out = computeGD(species = "Turnera_ulmifolia", gmatrix = 1, dmatrix = 1)
out$gmat

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=100)

save(modD, file="analyses/adj_Dmats/Turnera.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Turnera.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = out$G*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

vals = plot_GD(out$G, out$D, species="Turnera ulmifolia")

gdList[[length(gdList)+1]] = data.frame(sp = "Turnera_ulmifolia", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#### Fragaria ####

#NB: G-matrices are for females and hermaphrodites separately, D is for overall means
out = computeGD(species = "Fragaria_virginiana", gmatrix = 3, dmatrix = 1)
out$gmat

# Estimate error-corrected D matrix
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

modD = estimateDlog(means, eV, thin=100)

save(modD, file="analyses/adj_Dmats/Fragaria.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Fragaria.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = out$G*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

vals = plot_GD(out$G, out$D, species="Fragaria virginiana", plot=F)

gdList[[length(gdList)+1]] = data.frame(sp = "Fragaria_virginiana", g = out$gmat, traits = ncol(out$G), 
                                        emean = evolvabilityMeans(out$G)[1],
                                        emin = evolvabilityMeans(out$G)[2],
                                        emax = evolvabilityMeans(out$G)[3],
                                        cmean = evolvabilityMeans(out$G)[4],
                                        imean = evolvabilityMeans(out$G)[7],
                                        d = out$dmat, npops = out$nPop, 
                                        dmean = evolvabilityMeans(out$D)[1],
                                        betaG = vals[1,1], r2G = vals[2,1],
                                        betaD = vals[3,1], r2D = vals[4,1],
                                        theta = vals[5,1], row.names = NULL)

#### dgList output ####

save(gdList, file="gdList.RData")

load(file="gdList.RData")

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

#load(file="analyses/delph_Silene/gdDF.RData")
#gdList[[length(gdList)+1]] = gdDF

gdDat = rbind.fill(gdList)
dim(gdDat)
tail(gdDat)

for(i in 1: nrow(gdDat)){
  if(gdDat$theta[i]>90){
    gdDat$theta[i]=180-gdDat$theta[i]
  }
}

meanDat = ddply(gdDat, .(sp), summarize,
                emean = median(emean),
                cmean = median(cmean),
                dmean= median(dmean),
                betaG = median(betaG),
                r2G = median(r2G),
                betaD = median(betaD),
                r2D = median(r2D),
                theta = median(theta))
meanDat

min(c(range(gdDat$betaG), range(gdDat$betaD)))
max(c(range(gdDat$betaG), range(gdDat$betaD)))

x11(height=6, width=4.5)
mat = matrix(c(1,2,3,3,3,3),nrow=3, byrow=T)
layout(mat = mat)
par(mar=c(2,4,4,2))
hist(gdDat$r2G, xlab="", ylab="", main=expression(paste(R^2, " G eigenvectors")), las=1)
#legend("topleft", "R2 G", bty="n")
hist(gdDat$r2D, xlab="", ylab="", main=expression(paste(R^2, " D eigenvectors")), las=1)
#legend("topleft", "R2 D", bty="n")

par(mar=c(4,4,1,2))
plot(gdDat$betaG, gdDat$betaD, cex=gdDat$r2G*4, lwd=2, col="lightgrey",
     ylim=c(0,3), xlim=c(0,2), xlab="", ylab="", las=1)
points(meanDat$betaG, meanDat$betaD, pch=16)
points(median(meanDat$betaG), median(meanDat$betaD), pch=16, col="blue", cex=1.5)
abline(h=1, lty=2)
abline(v=1, lty=2)
mtext("Slope of G eigenvectors", 1, line=2.5)
mtext("Slope of D eigenvectors", 2, line=2.5)


par(mfrow=c(1,1))
plot(gdDat$r2G, gdDat$betaG, pch=16, las=1, ylim=c(0,3), ylab="Slope", xlab="R2 for G eigenvectors")
points(gdDat$r2G, gdDat$betaD)
legend("topleft", pch=c(1,16), legend=c("D", "G"))
abline(h=1)

plot(gdDat$r2G, gdDat$betaG, pch=16, las=1, ylim=c(0,3), ylab="Slope", xlab="R2")
points(gdDat$r2D, gdDat$betaD)
legend("topleft", pch=c(1,16), legend=c("D", "G"))
abline(h=1)

par(mfrow=c(2,2))
plot(gdDat$theta, gdDat$betaG, pch=16, ylim=c(0, 1.6), xlim=c(0,90),las=1)
coefs = summary(lm(betaG~theta, data=gdDat))$coef
x = 0:90
y = coefs[1,1] + coefs[2,1]*x
#lines(x, y)

plot(gdDat$theta, gdDat$r2G, ylim=c(0,1), xlim=c(0,90), pch=16, las=1)
coefs = summary(lm(r2G~theta, data=gdDat))$coef
x = 0:90
y = coefs[1,1] + coefs[2,1]*x
#lines(x, y)

plot(gdDat$theta, gdDat$betaD, pch=16, ylim=c(0, 3), xlim=c(0,90),las=1)
coefs = summary(lm(betaG~theta, data=gdDat))$coef
x = 0:90
y = coefs[1,1] + coefs[2,1]*x
#lines(x, y)

plot(gdDat$theta, gdDat$r2D, ylim=c(0,1), xlim=c(0,90), pch=16, las=1)
coefs = summary(lm(r2G~theta, data=gdDat))$coef
x = 0:90
y = coefs[1,1] + coefs[2,1]*x
#lines(x, y)



# All the studies ####
nG = NULL
nD = NULL
for(s in c(1:length(both_sp))){
  species = both_sp[s]  
  nG[s] = length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD[s] = length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
}
cbind(both_sp, nG, nD)

reslist = list()
for(s in 1:length(both_sp)){
  for(g in 1:nG[s]){
    for(d in 1:nD[s]){
      res = computeGD(species = both_sp[s], gmatrix = g, dmatrix = d)[c(1:2, 5:18)]
      reslist[length(reslist)+1] = as.data.frame(unlist(res)[c(1:16)])
    }
  }
}
length(reslist)
reslist[[20]]

resmat = as.data.frame(matrix(NA, ncol=16, nrow=length(reslist)))
for(i in 1:length(reslist)){
  resmat[i,1:2]=as.character(reslist[[i]][1:2])
  resmat[i,3:13]=round(as.numeric(as.character(reslist[[i]][3:13])),2)
  resmat[i,14:16]=as.numeric(as.character(reslist[[i]][14:16]))
}

colnames(resmat)=c("D", "G", "npop", "theta", "betaG", "r2G", "betaD", "r2D","betaT","r2T","i_mean","e_mean","d_mean","nBetaG","nBetaD","nBetaT")

for(i in 1:nrow(resmat)){
  if(resmat$theta[i]>=90){
    resmat$theta[i] = 180 - resmat$theta[i]
  }
}

names(resmat)
#resmat = resmat[,c(2,1,16,3, 12,13,4,5,6)]
head(resmat)
View(resmat)

plot(resmat$theta, resmat$betaG, pch=16)
plot(resmat$theta, resmat$r2G, pch=16)

resmat = resmat[resmat$npop>2,]
resmat = resmat[resmat$nBetaG>2,]
resmat = resmat[resmat$nBetaD>2,]
resmat = resmat[resmat$nBetaT>2,]

mean(resmat$betaG)
hist(resmat$r2G)
hist(resmat$betaG)

plot(resmat$betaT, resmat$betaG)
plot(resmat$r2T, resmat$r2G)

plot(resmat$theta, resmat$betaG)

plot(resmat$npop, resmat$betaG)
plot(resmat$npop, resmat$r2G, cex=resmat$nBetaG*.5)

plot(resmat$i_mean, resmat$r2G, cex=resmat$nBetaG*.5, xlim=c(0,1))

#### Error-corrected D original code ####

drop = which(is.na(rowSums(means[-1])))
if(length(drop)>0){
  means=means[-drop]
  eV=ev[-drop]
}

#Set prior
n = ncol(means)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1))

means[,1:ncol(means)] = apply(means[,1:ncol(means)], 2, function(x) x*10)
mev = melt(eV)$value*100
data = means
vars = paste0(colnames(means), collapse=", ")
vars
vars = paste0("c(",vars,") ~-1+trait")

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

mod<-MCMCglmm(as.formula(noquote(vars)),
              rcov = ~us(trait):units,
              mev = mev,
              data = means, 
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

modD = matrix(apply(mod$VCV, 2, median)[2:(1+(ncol(means))^2)], nrow=n)
colnames(modD) = rownames(modD) = colnames(means)
modD = meanStdG(modD, colMeans(means))
round(modD, 3)

