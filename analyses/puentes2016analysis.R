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

#####################################
#### - Estimating the D matrix - ####
#####################################
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
thin = 50
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

save(mod, file="./analyses/puentes2016/Dmat75k.RData")

#### G matrix for STUC population #### 
stuc = subset(Gdat, pop=="STUC")
stuc = droplevels(stuc)
stuc$cID=1:nrow(stuc)
stuc$animal=paste(as.character(stuc$sire), as.character(stuc$dam), stuc$cID, sep="_")

#Mean-scale and multiply by 10
stuc[,c(6,7,9,14)]=apply(stuc[,c(6,7,9,14)], 2, function(x) 10*x/mean(x, na.rm=T))

ped=subset(stuc, select = c("animal", "dam", "sire"))
parentped <- cbind(animal=unique(c(as.character(stuc$sire), as.character(stuc$dam))), dam=NA, sire=NA)
ped<-rbind(parentped, ped)
#ped <- MCMCglmm:: prunePed(ped, stuc$animal) 
head(ped)

invA = inverseA(ped)$Ainv

n = 4
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 10
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(petal.width.mm, petal.length.mm, flowers, rosette.size.cm) ~ -1+trait,
                random = ~us(trait):animal,
                rcov = ~us(trait):units,
                data = stuc, ginverse = list(animal = invA),
                family = rep("gaussian", n), prior = prior, 
                nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="Z:/data/Evoldiv/analyses/puentes2016/Gmat_STUC.RData")
#150k in 23.3 min

#### G matrix for each population 
a = Sys.time()

for(population in levels(Gdat$pop)){

reddat = subset(Gdat, pop==population)
reddat = droplevels(reddat)
reddat$cID = 1:nrow(reddat)
reddat$animal = paste(as.character(reddat$sire), as.character(reddat$dam), reddat$cID, sep="_")

#Mean-scale and multiply by 10
reddat[,c(6,7,9,14)] = apply(reddat[,c(6,7,9,14)], 2, function(x) 10*x/mean(x, na.rm=T))

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
thin = 10
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a = Sys.time()
mod <- MCMCglmm(c(petal.width.mm, petal.length.mm, flowers, rosette.size.cm) ~ -1+trait,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = reddat, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

filename = paste0("C:/data/Evoldiv/analyses/puentes2016/Gmat_", population,".RData")
save(mod, file=filename)
}
Sys.time()-a

#####################################
######## - Analyse G vs. D - ########
#####################################

#The D matrix
dat = Gdat

load(file="./analyses/puentes2016/Dmat75k.RData")
n=4
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
dmat

plot(mod$VCV[,3])

levels(Gdat$pop)
load(file="analyses/puentes2016/Gmat_SPIT.RData")
load(file="analyses/puentes2016/Gmat_STUC.RData")
load(file="analyses/puentes2016/Gmat_STUS.RData")
load(file="analyses/puentes2016/Gmat_VIS.RData")

plot(mod$VCV[,12])
summary(mod$VCV)

gmat = matrix(apply(mod$VCV, 2, median)[1:16], nrow=4)
colnames(gmat) = rownames(gmat) = c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
gmat

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)

isSymmetric.matrix(dmat)
#solve(dmat)
#dmat[lower.tri(dmat)] <- t(dmat)[lower.tri(dmat)]
gmat = round(gmat, 10)
dmat = round(dmat, 10)

#Compute eigenvectors etc.
g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
#var_d_g = diag(t(g_ev) %*% dmat %*% g_ev)

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e
#var_d_d = diag(t(d_ev) %*% dmat %*% d_ev)

#Compute summary stats
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

x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Arabidopsis lyrata: VIS", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#Plot
plot(var_g_g, var_d_g, xlim=c(-2,15), ylim=c(-2,15),
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Arabidopsis lyrata: G1", las=1)
points(var_g_d, var_d_d, pch=16)
points(diag(gmat), diag(dmat), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))


#Bending the D matrix
eigen(dmat)
lambda=diag(eigen(dmat)$values)
lambda
lambda[4,4]=0.000001
#lambda[4,4]=0.01

d2=d_ev%*%lambda%*%t(d_ev)
eigen(d2)

dmat
d2

dmat=d2

