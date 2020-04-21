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

#### - Estimating the P mean matrix - ####
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

# Check convergence
summary(mod$VCV)
plot(mod$VCV[,1])

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
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat) = rownames(gmat) = colnames(dat)[4:8]
gmat

# The D matrix
load(file="./analyses/andersson_crepis/Dmat75k.RData")

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

source("code/plot_GD.R")
vals = plot_GD(gmat, dmat, MeanP, species="Crepis tectorum", plot="c")

gdDF = data.frame(species="Crepis_tectorum", g = "Crepis tectorum: Visby", ntraits = ncol(gmat), 
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Crepis tectorum: All", npops = 54, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaG = vals$res[3,3], r2G = vals$res[3,4],
                  betaD = vals$res[4,3], r2D = vals$res[4,4],
                  betaD_cond = vals$res[5,3], r2D_cond = vals$res[5,4],
                  betaP = vals$res[6,3], r2P = vals$res[6,4],
                  betaP_cond = vals$res[7,3], r2P_cond = vals$res[6,4],
                  r2All = vals$res[8,4],
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
outdat = computeDelta(G=gmat/100, means=means, z0=z0)

deltaDF = data.frame(sp="Crepis_tectorum", g="Crepis tectorum: Visby", traits=ncol(gmat), 
                            d="Crepis tectorum: All", pop=rownames(means), 
                            emean=evolvabilityMeans(gmat)[1],
                            emin=evolvabilityMeans(gmat)[2],
                            emax=evolvabilityMeans(gmat)[3],
                            cmean=evolvabilityMeans(gmat)[4],
                            div=outdat[,1], edelta=outdat[,2], cdelta=outdat[,3], 
                            theta=outdat[,4], row.names=NULL)
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

