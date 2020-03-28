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
  Plist[[i]] = cov(sub[,3:12])
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
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = colnames(dat)[3:12]
dmat = dmat*100
dmat

plot(mod$VCV[,1])

# The G matrix
Gdat = read.csv("data/walter/Data_Exp2_Gvariance.csv")
head(Gdat)

pops = unique(Gdat$Type)
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

load(file="analyses/walter2018/Gmat_Dune.RData")
load(file="analyses/walter2018/Gmat_Head.RData")
load(file="analyses/walter2018/Gmat_Table.RData")
load(file="analyses/walter2018/Gmat_Wood.RData")

#plot(mod$VCV[,21])
#summary(mod$VCV)

n=10
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat)=rownames(gmat)=c("Height", "MSL_W", "SB", "MSD", "Area", "P2A2", "Circularity", "Nindents.Peri", "IndentWidthMean", "IndentDepthMean")
gmat

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)

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
     main="Senecio pimpinellifolium", las=1)
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
