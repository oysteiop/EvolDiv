#############################################################
#### Senecio pinnatifolius, data from Walter et al. 2018 ####
#############################################################

rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)

# D-matrix summary stats
Ddat = read.csv("data/walter/Data_Exp1_Dmatrix.csv")
head(Ddat)

popmeans = apply(Ddat[,3:12], 2, function(x){tapply(x, Ddat$Population, mean, na.rm=T)})
vars = apply(Ddat[,3:12], 2, function(x){tapply(x, Ddat$Population, var, na.rm=T)})
ns = apply(Ddat[,3:12], 2, function(x){tapply(x, Ddat$Population, length)})

#### Estimating the mean P matrix ####
#Ddat = Ddat[Ddat$Ecotype=="Dune_N",]
Ddat
pops = unique(Ddat$Population)
length(pops)
Plist = list()
for(i in 1:length(pops)){
  sub = Ddat[Ddat$Population==pops[i],]
  Plist[[i]] = meanStdG(cov(sub[,3:12]), colMeans(sub[,3:12]))
  
}

MeanP = apply(simplify2array(Plist), 1:2, mean)

#### Estimating the D matrix ####
dat = read.csv("data/walter/Data_Exp1_Dmatrix.csv")
head(dat)

dat[,c(3:12)] = apply(dat[,c(3:12)], 2, function(x) log(x))
head(dat)

#Ten traits
n = 10
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*5
prior <- list(R=list(V=diag(n), nu=n+0.002-1), 
              G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a = Sys.time()
mod <- MCMCglmm(c(VegHeight, MSL_W, SB, MSD, Area, P2A2, Circularity, Nindents.Peri, IndentWidth, IndentDepth) ~ -1+trait,
                random = ~us(trait):Population,
                rcov = ~us(trait):units,
                data = dat,
                family = rep("gaussian", n), prior = prior, 
                nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/walter2018/Dmat150k.RData")

#G-matrix summary stats
Gdat = read.csv("data/walter/Data_Exp2_Gvariance.csv")
head(Gdat)
popmeans=apply(Gdat[,8:17], 2, function(x){tapply(x, Gdat$Type,mean,na.rm=T)})
round(t(popmeans),3)

vars=apply(Gdat[,8:17], 2, function(x){tapply(x, Gdat$Type,var,na.rm=T)})
round(t(vars),4)

#### Estimate G matrices #### 
levels(Gdat$Type)

a=Sys.time()
for(type in levels(Gdat$Type)){
  
reddat = subset(Gdat, Type==type)
reddat = droplevels(reddat)
reddat$cID=1:nrow(reddat)
reddat$animal=paste(as.character(reddat$sire), as.character(reddat$dam), reddat$cID, sep="_")

#Mean-scale and multiply by 10
names(reddat)
reddat[,c(8:17)]=apply(reddat[,c(8:17)], 2, function(x) 10*x/mean(x, na.rm=T))

ped=subset(reddat, select = c("animal", "dam", "sire"))
parentped <- cbind(animal=unique(c(as.character(reddat$sire), as.character(reddat$dam))), dam=NA, sire=NA)
ped<-rbind(parentped, ped)
#ped <- MCMCglmm:: prunePed(ped, reddat$animal) 
#head(ped)

invA = inverseA(ped)$Ainv

n = 10
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

mod<-MCMCglmm(c(Height, MSL_W, SB, MSD, Area, P2A2, Circularity, Nindents.Peri, IndentWidthMean, IndentDepthMean) ~ -1+trait*Block,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = reddat, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

filename=paste0("analyses/walter2018/Gmat_", type,".RData")
save(mod, file=filename)
}
Sys.time()-a

######## Analyse G vs. D ########

# The D matrix
dat = read.csv("data/walter/Data_Exp1_Dmatrix.csv")

load(file="./analyses/walter2018/Dmat150k.RData")
n=10
dpost = mod$VCV[,1:(n*n)]*100
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)*100
colnames(dmat) = rownames(dmat) = colnames(dat)[3:12]
dmat

#plot(mod$VCV[,1])

# The G matrix
Gdat = read.csv("data/walter/Data_Exp2_Gvariance.csv")
head(Gdat)

pops = unique(Gdat$Type)
nFam = tapply(Gdat$sire, Gdat$Type, function(x) length(unique(x)))

glist = list()
n=10
for(i in 1:length(pops)){
  pop = pops[i]
  filename = paste0("analyses/walter2018/Gmat_", pop,".RData")
  load(filename)
  gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
  colnames(gmat)=rownames(gmat)=c("Height", "MSL_W", "SB", "MSD", "Area", "P2A2", "Circularity", "Nindents.Peri", "IndentWidthMean", "IndentDepthMean")
  glist[[i]] = gmat
}

gmat = apply(simplify2array(glist), 1:2, mean)

# Individual Gmatrices
pops=c("Dune", "Head", "Table", "Wood")
p=4
for(p in 1:length(pops)){
  pop = pops[p]
  filename=paste0("analyses/walter2018/Gmat_", pop, ".RData")
  load(filename)
#plot(mod$VCV[,21])
#summary(mod$VCV)

n=10
gpost = mod$VCV[,1:(n*n)]
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat)=rownames(gmat)=c("Height", "MSL_W", "SB", "MSD", "Area", "P2A2", "Circularity", "Nindents.Peri", "IndentWidthMean", "IndentDepthMean")
gmat

#evolvabilityMeans(gmat)
#evolvabilityMeans(dmat)

source("code/computeGD.R")
source("code/alignMat.R")
vals = computeGD(gmat, dpost, MeanP, species="Senecio pinnatifolius",
                 linearonly=F, SE=T, fixD=F, nSample=1000, nPop=16, nFam=nFam[p], plot=F)

#Uncertainty over the posterior
out = list()
for(i in 1:1000){
  sgmat = matrix(gpost[i,], nrow=n)
  sdmat = matrix(dpost[i,], nrow=n)
  #sdmat = dmat
  out[[i]] = computeGD(sgmat, sdmat, MeanP, species="")   
}

slopes = lapply(out, function(x) x$res$slope)
slopemean = apply(simplify2array(slopes), 1, median)
slopeSE = apply(simplify2array(slopes), 1, sd)

#vals$res$slope_MC = slopemean
#vals$res$SE = slopeSE

vals

name = paste0("Senecio pinnatifolius: ", pop)
gdDF = data.frame(species="Senecio_pinnatifolius", g = name, ntraits = ncol(gmat), 
                  nPop = 16, nFam = nFam[p],
                  dims = "area+coun+line+rati",
                  ndims = 4,
                  traitgroups = "veg",
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Senecio pinnatifolius: Walter et al 2018 greenhouse", nPop = 16, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaT = vals$res[1,3], betaT_SE = vals$res[1,5], r2T = vals$res[1,6],
                  betaT_cond = vals$res[2,3], betaT_cond_SE = vals$res[2,5], r2T_cond = vals$res[2,6],
                  betaG = vals$res[3,3], betaG_SE = vals$res[3,5], r2G = vals$res[3,6],
                  betaD = vals$res[4,3], betaD_SE = vals$res[4,5], r2D = vals$res[4,6],
                  betaD_cond = vals$res[5,3], betaD_cond_SE = vals$res[5,5], r2D_cond = vals$res[5,6],
                  betaP = vals$res[6,3], betaP_SE = vals$res[6,5], r2P = vals$res[6,6],
                  betaP_cond = vals$res[7,3], betaP_cond_SE = vals$res[7,5], r2P_cond = vals$res[7,6],
                  r2All = vals$res[8,6],
                  theta = vals$theta, row.names = NULL)
head(gdDF)

filename = paste0("analyses/walter2018/gdDF_", pop, ".RData")
save(gdDF, file=filename)
}


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
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Senecio pinnatifolius", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

# Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

# Ellipse plot
source("GDellipse.R")
GDellipse(dmat, gmat, xlim=c(-7,7), ylim=c(-7, 7), main="Senecio: Dune")


#### Divergence vectors ####
Gdat = read.csv("data/walter/Data_Exp2_Gvariance.csv")

pops = c("Dune", "Head", "Table", "Wood")
nFam = tapply(Gdat$sire, Gdat$Type, function(x) length(unique(x)))

p=1
for(p in 1:length(pops)){
  pop = pops[p]
  
# The G matrix
filename = paste0("analyses/walter2018/Gmat_", pop,".RData")
load(filename)

n=10
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat)=rownames(gmat)=c("Height", "MSL_W", "SB", "MSD", "Area", "P2A2", "Circularity", "Nindents.Peri", "IndentWidthMean", "IndentDepthMean")
signif(gmat, 2)

#Pop means
Ddat = read.csv("data/walter/Data_Exp1_Dmatrix.csv")
popmeans = apply(Ddat[,3:12], 2, function(x){tapply(x, Ddat$Population, mean, na.rm=T)})
means = popmeans
ecotype = factor(substr(rownames(popmeans), 1, 1))

Gdat = read.csv("data/walter/Data_Exp2_Gvariance.csv")
head(Gdat)
z0 = apply(Gdat[,8:17], 2, function(x){tapply(x, Gdat$Type, mean, na.rm=T)})
z0 = z0[p,] #Choose reference

source("code/computeDelta.R")
outdat = computeDelta3(gmat/100, means, z0, SE=T, nFam=nFam[p], nSample=1000)

name = paste0("Senecio pinnatifolius: ", pop)
deltaDF = data.frame(species="Senecio_pinnatifolius", g=name, ntraits=ncol(gmat), 
                     d= "Senecio pinnatifolius: Walter et al 2018 greenhouse", pop=rownames(means), 
                     emean=outdat$emean,
                     emin=outdat$emin,
                     emax=outdat$emax,
                     cmean=outdat$cmean,
                     div=outdat$div,
                     edelta=outdat$edelta, edeltaSE=outdat$edeltaSE,
                     cdelta=outdat$cdelta, cdeltaSE=outdat$cdeltaSE,
                     theta=outdat$theta, row.names=NULL)
head(deltaDF)

filename=paste0("analyses/walter2018/deltaDF_", pop,".RData")

save(deltaDF, file=filename)
}

#save(deltaDF, file="analyses/walter2018/deltaDF_Dune.RData")
#save(deltaDF, file="analyses/walter2018/deltaDF_Head.RData")
#save(deltaDF, file="analyses/walter2018/deltaDF_Table.RData")
#save(deltaDF, file="analyses/walter2018/deltaDF_Wood.RData")


#x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0,12), las=1,
     xlab="", ylab="",
     main="Senecio pinnatifolius", col=c("black", "red", "blue", "green")[as.numeric(ecotype)])
mtext("Divergence from focal population(%)", 1, line=2.5)
mtext("Evolvability(%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 127
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

