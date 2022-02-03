#####################################################
#### Crepis tectorum, data from Stefan Andersson ####
#####################################################

rm(list=ls())
library(evolvability)
library(plyr)
library(reshape2)
library(MCMCglmm)
library(lme4)
library(kinship2)

list.files(path="./data/andersson")
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
head(dat)
str(dat)

#### - Estimating the mean P matrix - ####
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
dat = dat[dat$TYPE=="POP",]
head(dat)

pops = unique(dat$IDENTITY)
length(pops)

Plist = list()
for(i in 1:length(pops)){
  sub = dat[dat$IDENTITY==pops[i],]
  Plist[[i]] = meanStdG(cov(sub[,4:8]), colMeans(sub[,4:8]))
}

MeanP = apply(simplify2array(Plist), 1:2, mean)
MeanP

#### - Estimating the G matrix - ####
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
dat = dat[dat$TYPE=="OUTX",]
dat$animal = 1:nrow(dat)
head(dat)

# Build the pedigree
pedigree = data.frame(as.character(dat$animal))
pedigree$dam = paste0(dat$IDNO, "d")
pedigree$sire = paste0(dat$IDNO, "s")
names(pedigree) = c("animal","dam","sire")

parentped = data.frame(animal=c(unique(pedigree$dam), unique(pedigree$sire)))
parentped$dam = rep(NA, nrow(parentped))
parentped$sire = rep(NA, nrow(parentped))
names(parentped) = c("animal","dam","sire")

pedigree = rbind(parentped, pedigree)
head(pedigree)

# Mean-scale and multiply by 10
dat[,c(4:8)] = apply(dat[,c(4:8)], 2, function(x) 10*x/mean(x, na.rm=T))
head(dat)

# Run the MCMCglmm analysis
#Aped <- 2 * kinship2::kinship(pedigree[, 1], pedigree[,2], pedigree[, 3])

invA = inverseA(pedigree)$Ainv

# Five traits
n = 5
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior <- list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(LEN,TIP,MAX,MIN,TEETH) ~ -1+trait,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = dat, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/andersson_crepis/Gmat75k.RData")

# Simpler model
n = 5
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior <- list(R=list(V=diag(n), nu=n+0.002-1), 
              G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 25
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(LEN,TIP,MAX,MIN,TEETH) ~ -1+trait,
              random = ~us(trait):IDENTITY,
              rcov = ~us(trait):units,
              data = dat,
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/andersson_crepis/Gmat75k_FS.RData")

# Check convergence
summary(mod$VCV)
plot(mod$VCV[,1])

# HMSC
library(Hmsc)
names(dat)
Y = data.frame(dat[,4:8])
dfPi = data.frame(dat[,2])
names(dfPi)="IDENTITY"
dfPi$IDENTITY=factor(dfPi$IDENTITY)

rL1 = HmscRandomLevel(units=unique(dfPi$IDENTITY))

rL1$nfMin=5
rL1$nfMax=5

m = Hmsc(Y = as.matrix(Y), XData=data.frame(rep(1, nrow(Y))), XFormula = ~1,  dist = "normal", 
         studyDesign = dfPi, ranLevels=list(IDENTITY=rL1))

# RUN MCMC
samples = 1000
thin = 1
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 1

a1 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains, updater=list(GammaEta=FALSE))
b1 = Sys.time()

#2*300k in 7.4 days
b1-a1 #15k in 22 min

save(m, file="./analyses/andersson_crepis/Gmat15k_Hmsc.RData")

# Check convergence
post = convertToCodaObject(m)
str(post$Lambda)
plot(post$Lambda[[1]][,1:2])

Omegapost = getPostEstimate(m, "Omega")
gest = Omegapost$mean*2
gest

#### - Estimating the D matrix - ####
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
dat = dat[dat$TYPE=="POP",]
head(dat)

dat[,c(4:8)] = apply(dat[,c(4:8)], 2, function(x) log(x))
head(dat)

# Five traits
n = 5
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*5
prior <- list(R=list(V=diag(n), nu=n+0.002-1), 
              G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a = Sys.time()
mod <- MCMCglmm(c(LEN,TIP,MAX,MIN,TEETH) ~ -1+trait,
              random = ~us(trait):IDENTITY,
              rcov = ~us(trait):units,
              data = dat,
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/andersson_crepis/Dmat75k.RData")

# Check convergence
summary(mod$VCV)
plot(mod$VCV[,19])

#### Divergence analysis ####
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")

# The G matrix
load(file="./analyses/andersson_crepis/Gmat75k.RData")
n = 5
gpost = mod$VCV[, 1:(n*n)]
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat) = rownames(gmat) = colnames(dat)[4:8]
gmat

# The D matrix
load(file="./analyses/andersson_crepis/Dmat75k.RData")

dpost = mod$VCV[,1:(n*n)]*100
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = colnames(dat)[4:8]
dmat = dmat*100
dmat

# EvolvabilityMeans
evolvabilityMeans(gmat)
#evolvabilityMeansMCMC(mod$VCV[,1:(n*n)])
evolvabilityMeans(dmat)
signif(cov2cor(gmat), 2)
signif(cov2cor(dmat), 2)

source("code/computeGD.R")
source("code/alignMat.R")
vals = computeGD(gmat, dpost, MeanP, species="Crepis tectorum", plot=F)

# Uncertainty over the posterior
out = list()
for(i in 1:10){
  sgmat = matrix(gpost[i,], nrow=n)
  sdmat = matrix(dpost[i,], nrow=n)
  #sdmat = dmat
  out[[i]] = computeGD(sgmat, sdmat, MeanP, species="")   
}

slopes = lapply(out, function(x) x$res$slope)
slopemean = apply(simplify2array(slopes), 1, median)
slopeSE = apply(simplify2array(slopes), 1, sd)

vals$res$slope_MC = slopemean
vals$res$SE = slopeSE

vals

gdDF = data.frame(species="Crepis_tectorum", g = "Crepis tectorum: Vickleby", ntraits = ncol(gmat), 
                  dims = "line+coun",
                  ndims = 2,
                  traitgroups = "veg",
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Crepis tectorum: Andersson Crepis greenhouse", nPop = 54, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaT = vals$res[1,3], betaT_SE = vals$res[1,5], r2T = vals$res[1,6],
                  betaT_cond = vals$res[2,3], r2T_cond = vals$res[2,6],
                  betaG = vals$res[3,3], betaG_SE = vals$res[3,5], r2G = vals$res[3,6],
                  betaD = vals$res[4,3], betaD_SE = vals$res[4,5], r2D = vals$res[4,6],
                  betaD_cond = vals$res[5,3], r2D_cond = vals$res[5,6],
                  betaP = vals$res[6,3], r2P = vals$res[6,6],
                  betaP_cond = vals$res[7,3], r2P_cond = vals$res[7,6],
                  r2All = vals$res[8,6],
                  theta = vals$theta, row.names = NULL)
head(gdDF)

save(gdDF, file="analyses/andersson_crepis/gdDF.RData")

# Compute eigenvectors etc.
g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e

p_ev = eigen(MeanP)$vectors
var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e

# Compute summary stats
mt = lm(log(diag(dmat))~log(diag(gmat)))
beta_t = summary(mt)$coef[2,1]
beta_t
r2_t = summary(mt)$r.squared
r2_t

mg = lm(log(var_d_g)~log(var_g_g))
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d))
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

mp = lm(log(var_d_p)~log(var_g_p))
beta_p = summary(mp)$coef[2,1]
beta_p
r2_p = summary(mp)$r.squared
r2_p

# Plot
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Crepis tectorum", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")

legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

# Plot with modifed axes

cairo_pdf("multiEvolDivCases.pdf", height=8, width=8, fam="Times")
par(mfrow=c(2,2), mar=c(4,4,2,2))

#x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), pch=16, 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="", 
     ylab="",
     xaxt="n",
     yaxt="n",
     main="", las=1)
mtext("Proportional divergence", 2, line=3.0, cex=1)
mtext(expression(paste(italic(Crepis), " ", italic(tectorum))), line=0.5, cex=.9)

points(log10(var_g_d), log10(var_d_d), pch=1)
points(log10(var_g_p), log10(var_d_p), pch=16, col="firebrick")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue3")

mean1 = mean(log10(c(diag(dmat), var_d_g, var_d_d)))
mean2 = mean(log10(c(diag(gmat), var_g_g, var_g_d)))
segments(x0=mean2-10, y0=mean1-10, x1=mean2+10, y1=mean1+10)

legend("bottomright", legend=c(paste0("Original traits (", round(100*r2_t, 1),"%)"),
                               paste0("G-directions (", round(100*r2_g, 1),"%)"),
                               paste0("D-directions (", round(100*r2_d, 1),"%)"),
                               paste0("P-directions (", round(100*r2_p, 1),"%)")),
       pch=c(16, 16, 1, 16), col=c("blue3", "black", "black", "firebrick"))

axis(1, at=c(-.5, 0, .5, 1, 1.5), signif(10^c(-.5, 0, .5, 1, 1.5),1))

x3at = seq(-.5, 1.5, .5)
x3 = exp(sqrt(((10^x3at)/100)*(2/pi)))
axis(2, at=x3at, signif(x3, 3), las=1)

#xt3 = c(1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
#x3at = log10(100*log(xt3)^2/(2/pi))
#axis(2, at=x3at, signif(xt3, 4), las=1)

####

cvals=NULL
for(i in 1:5){
b=rep(0,5)
b[i]=1
cvals[i] = evolvabilityBeta(gmat, b)$c
}
points(log10(cvals), log10(diag(dmat)), pch=16, col="blue")

# Angles
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)





#### Divergence vectors ####
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")

# The G matrix
load(file="./analyses/andersson_crepis/Gmat75k.RData")
n = 5
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat) = rownames(gmat) = colnames(dat)[4:8]
round(gmat, 2)

# Pop means
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
dat = dat[dat$TYPE=="POP",]
dat$IDENTITY = factor(dat$IDENTITY)
means = apply(dat[,4:8], 2, function(x) tapply(x, dat$IDENTITY, mean, na.rm=T))

head(means)
dim(means)

dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
dat = dat[dat$TYPE=="OUTX",]
z0 = colMeans(dat[,4:8])

source("code/computeDelta.R")
outdat = computeDelta2(G=gmat/100, means=means, z0=z0)

deltaDF = data.frame(species="Crepis_tectorum", g="Crepis tectorum: Vickleby", ntraits=ncol(gmat), 
                            d="Crepis tectorum: Andersson Crepis greenhouse", pop=rownames(means), 
                            emean=outdat$emean,
                            emin=outdat$emin,
                            emax=outdat$emax,
                            cmean=outdat$cmean,
                            div=outdat$div, edelta=outdat$edelta, cdelta=outdat$cdelta,
                            theta=outdat$theta, row.names=NULL)

head(deltaDF)

save(deltaDF, file="analyses/andersson_crepis/deltaDF.RData")

#x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0,50), las=1,
     xlab="", ylab="",
     main="Crepis tectorum")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 87
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

