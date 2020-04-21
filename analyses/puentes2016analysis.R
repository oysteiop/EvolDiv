#############################################################
#### Arabidopsis lyrata, data from Puentes et al. 2016 - ####
#############################################################

rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)

Gdat = read.csv("data/puentes/puentes_etal_phen_multivariate_G_matrices.csv")
Gdat$rosette.size.cm = sqrt(Gdat$rosette.size.cm2)
head(Gdat)

popmeans = apply(Gdat[,5:13], 2, function(x) tapply(x, Gdat$pop, mean, na.rm = T))
popmeans = popmeans[,c(2,3,5,9)]
popmeans

vars = apply(Gdat[,5:13], 2, function(x) tapply(x, Gdat$pop, var, na.rm=T))
vars = vars[,c(2,3,5,9)]
vars

#### - Estimating the average P matrix - ####
dat = Gdat
head(dat)

pops = unique(dat$pop)
length(pops)
Plist = list()
for(i in 1:length(pops)){
        sub = dat[dat$pop==pops[i],]
        Plist[[i]] = meanStdG(cov(sub[,c(6,7,9,13)]), colMeans(sub[,c(6,7,9,13)]))
}

MeanP = apply(simplify2array(Plist), 1:2, mean)
MeanP

#### - Estimating the D matrix - ####
dat = Gdat
head(dat)
dat$fam = paste(dat$dam, dat$sire, sep="_")

#Log-transform and multiple by 10
dat[,c(6,7,9,13)] = apply(dat[,c(6,7,9,13)], 2, function(x) log(x)*10)
head(dat)

#Four traits
n = 4
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*5
prior <- list(R=list(V=diag(n), nu=n+0.002-1), 
              G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                     G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a = Sys.time()
mod <- MCMCglmm(c(petal.width.mm, petal.length.mm, flowers, rosette.size.cm) ~ -1+trait,
                random = ~us(trait):pop + us(trait):fam,
                rcov = ~us(trait):units,
                data = dat,
                family = rep("gaussian", n), prior = prior, 
                nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/puentes2016/Dmat150k.RData")

#### Estimate G matrices for each population #### 
a = Sys.time()

for(population in levels(Gdat$pop)){

reddat = subset(Gdat, pop==population)
reddat = droplevels(reddat)
reddat$cID = 1:nrow(reddat)
reddat$animal = paste(as.character(reddat$sire), as.character(reddat$dam), reddat$cID, sep="_")

# Mean-scale and multiply by 10
reddat[,c(6,7,9,13)] = apply(reddat[,c(6,7,9,13)], 2, function(x) 10*x/mean(x, na.rm=T))

ped = subset(reddat, select = c("animal", "dam", "sire"))
parentped = cbind(animal=unique(c(as.character(reddat$sire), as.character(reddat$dam))), dam=NA, sire=NA)
ped = rbind(parentped, ped)

invA = inverseA(ped)$Ainv

n = 4
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

mod <- MCMCglmm(c(petal.width.mm, petal.length.mm, flowers, rosette.size.cm) ~ -1+trait,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = reddat, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

filename = paste0("analyses/puentes2016/Gmat_", population,".RData")
save(mod, file=filename)
}
Sys.time()-a

######## Analyse G vs. D ########

# The D matrix
dat = Gdat

load(file="./analyses/puentes2016/Dmat150k.RData")
n=4
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
dmat

#plot(mod$VCV[,3])
pops = c("SPIT", "STUC", "STUS", "VIS")
glist = list()
for(i in 1:length(pops)){
pop = pops[i]
filename = paste0("analyses/puentes2016/Gmat_", pop,".RData")
load(filename)
gmat = matrix(apply(mod$VCV, 2, median)[1:16], nrow=4)
colnames(gmat) = rownames(gmat) = c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
glist[[i]] = gmat
}

gmat = apply(simplify2array(glist), 1:2, mean)

levels(Gdat$pop)
load(file="analyses/puentes2016/Gmat_SPIT.RData")
load(file="analyses/puentes2016/Gmat_STUC.RData")
load(file="analyses/puentes2016/Gmat_STUS.RData")
load(file="analyses/puentes2016/Gmat_VIS.RData")

#plot(mod$VCV[,1])
#summary(mod$VCV)

gmat = matrix(apply(mod$VCV, 2, median)[1:16], nrow=4)
colnames(gmat) = rownames(gmat) = c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
gmat

source("code/plot_GD.R")
vals = plot_GD(gmat, dmat, MeanP, species="Arabidopsis lyrata", plot=F)

gdDF = data.frame(species="Arabidopsis_lyrata", g = "Arabidopsis lyrata: VIS", ntraits = ncol(gmat), 
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Arabidopsis lyrata: All", npops = 4, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaG = vals$res[3,3], r2G = vals$res[3,4],
                  betaD = vals$res[4,3], r2D = vals$res[4,4],
                  betaD_cond = vals$res[5,3], r2D_cond = vals$res[5,4],
                  betaP = vals$res[6,3], r2P = vals$res[6,4],
                  betaP_cond = vals$res[7,3], r2P_cond = vals$res[6,4],
                  r2All = vals$res[8,4],
                  theta = vals$theta, row.names = NULL)
head(gdDF)

save(gdDF, file="analyses/puentes2016/gdDF_SPIT.RData")
save(gdDF, file="analyses/puentes2016/gdDF_STUC.RData")
save(gdDF, file="analyses/puentes2016/gdDF_STUS.RData")
save(gdDF, file="analyses/puentes2016/gdDF_VIS.RData")



evolvabilityMeans(gmat)
evolvabilityMeans(dmat)
signif(cov2cor(gmat), 2)
signif(cov2cor(dmat), 2)

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

x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [x100])", 
     ylab="log10 (Divergence [%])", 
     main="Arabidopsis lyrata: VIS", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))


#### Divergence vectors ####

# The G matrix
load(file="analyses/puentes2016/Gmat_SPIT.RData")
load(file="analyses/puentes2016/Gmat_STUC.RData")
load(file="analyses/puentes2016/Gmat_STUS.RData")
load(file="analyses/puentes2016/Gmat_VIS.RData")

gmat = matrix(apply(mod$VCV, 2, median)[1:16], nrow=4)
colnames(gmat) = rownames(gmat) = c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
round(gmat, 2)

Gdat = read.csv("data/puentes/puentes_etal_phen_multivariate_G_matrices.csv")
Gdat$rosette.size.cm = sqrt(Gdat$rosette.size.cm2)
popmeans = apply(Gdat[,5:13], 2, function(x) tapply(x, Gdat$pop, mean, na.rm = T))
popmeans = popmeans[,c(2,3,5,9)]
z0 = popmeans[4,] #Select focal pop

# Pop means
means = popmeans[-4, ] #Select focal pop
head(means)
dim(means)

outdat = computeDelta(gmat/100, means, z0)

deltaDF = data.frame(sp="Arabidopsis_lyrata", g="Arabidopsis lyrata: VIS", traits=ncol(gmat), 
                     d="Arabidopsis lyrata: All", pop=rownames(means), 
                     emean=evolvabilityMeans(gmat)[1],
                     emin=evolvabilityMeans(gmat)[2],
                     emax=evolvabilityMeans(gmat)[3],
                     cmean=evolvabilityMeans(gmat)[4],
                     div=outdat[,1], edelta=outdat[,2], cdelta=outdat[,3], 
                     theta=outdat[,4], row.names=NULL)
head(deltaDF)

save(deltaDF, file="analyses/puentes2016/deltaDF_SPIT.RData")
save(deltaDF, file="analyses/puentes2016/deltaDF_STUC.RData")
save(deltaDF, file="analyses/puentes2016/deltaDF_STUS.RData")
save(deltaDF, file="analyses/puentes2016/deltaDF_VIS.RData")



x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0, 10), las=1,
     xlab="", ylab="",
     main="Arabidopsis lyrata: SPIT")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 33
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")
