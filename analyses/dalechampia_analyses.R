##############################################
#### Dalechampia 'scandensoides' analyses ####
##############################################

rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)

###### G matrix for the Tulum population (Mexico) ######
dat = read.csv("data/dalechampia/tulum/TulumDiallel.csv", header=T)
head(dat)

ped = read.csv("data/dalechampia/tulum/Pedigree.csv", header=T)
head(ped, 5)

names(dat)

# Define some traits
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$animal = dat$ind.

nms = c("BA", "GA", "GSD", "SW", "GAD", "ASD")

# G matrix
X = subset(dat, select=c(BA, GA, GSD, SW, GAD, ASD, animal, ind., Dam, meas.date, stage))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 100*x/mean(x, na.rm=T))
X = na.omit(X)

ped2 = MCMCglmm::prunePed(ped, X$animal)
invA = inverseA(ped2)$Ainv

Pres = var(X[,1:length(nms)], na.rm=T)-var(apply(X[,1:length(nms)], 2, function(x) tapply(x, X$animal, mean, na.rm=T)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G3=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(BA, GA, GSD, SW, GAD, ASD) ~ -1+trait:factor(stage),
              random= ~us(trait):animal + us(trait):ind. + us(trait):meas.date,
              rcov=~us(trait):units,
              data=X, ginverse=list(animal = invA),
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/dalechampia/Gmat_Tulum.RData")

###### G matrix for the Tovar population (Venezuela) ######
dat = read.csv("data/dalechampia/tovar/Merida diallel.csv", header=T)
head(dat)
names(dat)

# Define some traits
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$Gwtot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$animal = dat$ind

nms = c("BA", "GA", "GSD", "SW", "GAD")

# G matrix
X = subset(dat, select=c(BA, GA, GSD, SW, GAD, animal, ind, dam, measdate, male.fl))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 100*x/mean(x, na.rm=T))
X = na.omit(X)

# Constructing pedigree
dt2 = read.csv("data/dalechampia/tovar/mother identity.csv")

ped <- cbind(animal=unique(c(dat$sire, dat$dam)), dam=NA, sire=NA)
for(i in 1:nrow(ped)) ped[i, "dam"] <- as.character(dt2$X[which(dt2$Code == ped[i,"animal"])])

ped <- rbind(
        cbind(animal=unique(as.character(dt2$X)), dam=NA, sire=NA),
        ped,
        with(dat, cbind(animal=ind, dam=dam, sire=sire))
)

ped<-as.data.frame(ped)
ped2 <- MCMCglmm:: prunePed(ped, dat$ind) 
invA = inverseA(ped2)$Ainv

Pres = var(X[,1:length(nms)], na.rm=T)-var(apply(X[,1:length(nms)], 2, function(x) tapply(x, X$animal, mean, na.rm=T)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G3=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(BA, GA, GSD, SW, GAD) ~ -1+trait:factor(male.fl),
              random= ~us(trait):animal + us(trait):ind + us(trait):measdate,
              rcov=~us(trait):units,
              data=X, ginverse=list(animal = invA),
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/dalechampia/Gmat_Tovar.RData")

###### D and P matrix for Mexican populations of D. scandensoides ######
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)
head(dat)

# Define some traits
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$ind = paste(dat$pop, dat$plant, dat$seed, sep="_")
dat$ind = factor(dat$ind)
dat$population = factor(dat$population)

nms = c("BA", "GA", "GSD", "SW", "GAD", "ASD")

unique(dat$population)

pdat = dat[dat$population %in% c("BA", "C", "CC", "GS", "LM", "M", "PM"),]
pdat = na.omit(subset(pdat, select=c(BA, GA, GSD, SW, GAD, ASD, population)))
table(pdat$population)

# Mean P matrix
pops = unique(pdat$population)
length(pops)
Plist = list()
for(i in 1:length(pops)){
  sub = pdat[pdat$population==pops[i],]
  pm = cov(sub[,1:6])
  Plist[[i]] = meanStdG(pm, colMeans(sub[,1:6]))
}

MeanP = apply(simplify2array(Plist), 1:2, mean)
save(MeanP, file="analyses/dalechampia/MeanP_MX.RData")

# D matrix
dat = dat[dat$population %in% c("BA", "C", "CA", "CC", "CN", "CP", "GS", "LM", "M", "PM", "T", "Tulum51"),]
dat = dat[dat$ASD>0,]

X = subset(dat, select=c(BA, GA, GSD, SW, GAD, ASD, population, ind))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 10*log(x))
X = na.omit(X)
head(X)
str(X)

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(BA, GA, GSD, SW, GAD, ASD) ~ -1+trait,
              random= ~us(trait):population + us(trait):ind,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/dalechampia/Dmat_Mexico_75k.RData")

###### D matrix for Mexican populations of D. scandens ######
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)
head(dat)

# Define some traits
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$ind = paste(dat$pop, dat$plant, dat$seed, sep="_")
dat$ind = factor(dat$ind)
dat$population = factor(dat$population)

nms = c("BA", "GA", "GSD", "SW", "GAD")

unique(dat$population)

pdat = dat[dat$population %in% c("CH", "CO", "E", "HO", "P", "S", "V", "XU"),]
pdat = na.omit(subset(pdat, select=c(BA, GA, GSD, SW, GAD, population)))
table(factor(pdat$population))

# Mean P matrix
pops = unique(pdat$population)
length(pops)
Plist = list()
for(i in 1:length(pops)){
  sub = pdat[pdat$population==pops[i],]
  pm = cov(sub[,1:5])
  Plist[[i]] = meanStdG(pm, colMeans(sub[,1:5]))
}

MeanP = apply(simplify2array(Plist), 1:2, mean)
save(MeanP, file="analyses/dalechampia/MeanP_SG.RData")

# D matrix
unique(dat$population)
dat = dat[dat$population %in% c("CH", "CO", "D", "E", "HO", "MO", "P", "S", "V", "X", "XU", "Tovar70"),]

X = subset(dat, select=c(BA, GA, GSD, SW, GAD, population, ind))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 10*log(x))
X = na.omit(X)
head(X)
str(X)

X$population=factor(X$population)
write.csv2(melt(t(apply(X[,1:5], 2, function(x) tapply(x>0, X$population, sum, na.rm=T)))), "SGmeans.csv")


samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(BA, GA, GSD, SW, GAD) ~ -1+trait,
              random= ~us(trait):population + us(trait):ind,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/dalechampia/Dmat_Mexico_SG_75k.RData")

###### D and P matrix for Costa Rican populations ######
dat = read.table("data/dalechampia/populations/costarica_greenhouse.txt", header=T)
head(dat)

# Define some traits
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GW)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$ind = dat$ID
dat$ind = factor(dat$ind)
dat$population = factor(dat$Pop)

nms = c("BA", "GA", "GSD", "SW", "GAD", "ASD")

unique(dat$population)
length(unique(dat$population))
table(dat$population)

pdat = dat[dat$population %in% c("S1","S2","S6","S7","S8","S9","S11","S12","S20","S21","S22","S26"),]
pdat = na.omit(subset(pdat, select=c(BA, GA, GSD, SW, GAD, ASD, population)))
#pdat = pdat[pdat$ASD>0,]

# Mean P matrix
pops = unique(pdat$population)
length(pops)
Plist = list()
for(i in 1:length(pops)){
  sub = pdat[pdat$population==pops[i],]
  pm = cov(sub[,1:6])
  Plist[[i]] = meanStdG(pm, colMeans(sub[,1:6]))
}

MeanP = apply(simplify2array(Plist), 1:2, mean)
save(MeanP, file="analyses/dalechampia/MeanP.RData")

# D matrix
dat = dat[dat$population %in% c("S1","S2","S3","S6","S7","S8","S9","S11","S12","S13","S16","S18","S19","S20","S21","S22","S26"),]
dat = dat[dat$ASD>0,]

X = subset(dat, select=c(BA, GA, GSD, SW, GAD, ASD, population, ind))
X[,1:length(nms)] = apply(X[,1:length(nms)], 2, function(x) 10*log(x))
X = na.omit(X)
head(X)
str(X)

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(nms)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(BA, GA, GSD, SW, GAD, ASD) ~ -1+trait,
              random= ~us(trait):population + us(trait):ind,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/dalechampia/Dmat_CostaRica_75k.RData")

######## Analyse G vs. D: Mexico ########

#The D matrix
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)

load(file="analyses/dalechampia/Dmat_Mexico_75k.RData")
n=6
dpost = mod$VCV[,1:(n*n)]
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
dmat

#plot(mod$VCV[,1])

#The G matrix
dat = read.csv("data/dalechampia/tulum/TulumDiallel.csv", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tulum.RData")

n = 6
gpost = mod$VCV[,1:(n*n)]/100
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

#The P matrix
load(file="analyses/dalechampia/MeanP_MX.Rdata")
pmat = MeanP

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)
evolvabilityMeans(pmat*100)

source("code/computeGD.R")
source("code/alignMat.R")
vals = computeGD(gmat, dmat, pmat, species="Dalechampia scandens A", plot=F)

#Uncertainty over the posterior
out = list()
for(i in 1:100){
  sgmat = matrix(gpost[i,], nrow=n)
  sdmat = matrix(dpost[i,], nrow=n)
  #sdmat = dmat
  out[[i]] = computeGD(sgmat, sdmat, pmat, species="Dalechampia scandens A")   
}

slopes = lapply(out, function(x) x$res$slope)
slopemean = apply(simplify2array(slopes), 1, median)
slopeSE = apply(simplify2array(slopes), 1, sd)

vals$res$slope_MC = slopemean
vals$res$SE = slopeSE

vals

gdDF = data.frame(species="Dalechampia_scandens_A", g = "Dalechampia scandens A: Tulum", ntraits = ncol(gmat), 
                  dims = "line",
                  ndims = 1,
                  traitgroups = "flo",
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Dalechampia scandens A: Bolstad et al. Mexico greenhouse" , nPop = 12, 
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

save(gdDF, file="analyses/dalechampia/gdDF_Tulum_MX.RData")

# Compute eigenvectors etc.
g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e

p_ev = eigen(pmat)$vectors
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
     ylab="log10 (Divergence [%])", 
     main="Dalechampia: Tulum Mexico", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

######## Analyse G vs. D: Costa Rica ########

# The D matrix
dat = read.table("data/dalechampia/populations/costarica_greenhouse.txt", header=T)
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GW)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)

load(file="analyses/dalechampia/Dmat_CostaRica_75k.RData")
n=6
dpost = mod$VCV[,1:(n*n)]
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
dmat

#plot(mod$VCV[,1])

# The G matrix
dat = read.csv("data/dalechampia/tulum/TulumDiallel.csv", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tulum.RData")

n = 6
gpost = mod$VCV[,1:(n*n)]
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

#The mean P matrix
load(file="analyses/dalechampia/MeanP.RData")
pmat = MeanP

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)
evolvabilityMeans(pmat)
eigen(gmat)
eigen(dmat)
eigen(pmat)

source("code/computeGD.R")
source("code/alignMat.R")
vals = computeGD(gmat, dmat, pmat, species="Dalechampia scandens", ymin=-1, plot=F)

# Uncertainty over the posterior
out = list()
for(i in 1:10){
  sgmat = matrix(gpost[i,], nrow=n)
  sdmat = matrix(dpost[i,], nrow=n)
  #sdmat = dmat
  out[[i]] = computeGD(sgmat, sdmat, pmat, species="Dalechampia scandens A")   
}

slopes = lapply(out, function(x) x$res$slope)
slopemean = apply(simplify2array(slopes), 1, median)
slopeSE = apply(simplify2array(slopes), 1, sd)

vals$res$slope_MC = slopemean
vals$res$SE = slopeSE

vals

gdDF = data.frame(species="Dalechampia_scandens_A", g = "Dalechampia scandens A: Tulum", ntraits = ncol(gmat), 
                  dims = "line",
                  ndims = 1,
                  traitgroups = "flo",
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Dalechampia scandens A: Opedal et al. Costa Rica greenhouse", nPop = 17, 
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

save(gdDF, file="analyses/dalechampia/gdDF_Tulum_CR.RData")

# Compute eigenvectors etc.
g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_g_g_c = evolvabilityBeta(gmat, Beta = g_ev)$c
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_g_d_c = evolvabilityBeta(gmat, Beta = d_ev)$c
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e

p_ev = eigen(pmat)$vectors
var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
var_g_p_c = evolvabilityBeta(gmat, Beta = p_ev)$c
var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e

# Compute summary stats
cvals=NULL
for(i in 1:ncol(gmat)){
  b=rep(0,ncol(gmat))
  b[i]=1
  cvals[i] = evolvabilityBeta(gmat, b)$c
}

mt = lm(log(diag(dmat))~log(diag(gmat)))
beta_t = summary(mt)$coef[2,1]
beta_t
r2_t = summary(mt)$r.squared
r2_t

mtc = lm(log(diag(dmat))~log(cvals))
beta_tc = summary(mtc)$coef[2,1]
beta_tc
r2_tc = summary(mtc)$r.squared
r2_tc

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

mdc = lm(log(var_d_d)~log(var_g_d_c))
beta_dc = summary(mdc)$coef[2,1]
beta_dc
r2_dc = summary(mdc)$r.squared
r2_dc

mp = lm(log(var_d_p)~log(var_g_p))
beta_p = summary(mp)$coef[2,1]
beta_p
r2_p = summary(mp)$r.squared
r2_p

mpc = lm(log(var_d_p)~log(var_g_p_c))
beta_pc = summary(mpc)$coef[2,1]
beta_pc
r2_pc = summary(mpc)$r.squared
r2_pc

res = data.frame(traits = c("Original", "Original","G eigenvectors", "D eigenvectors", "D eigenvectors", "P eigenvectors", "P eigenvectors"),
                 evol = c("e","c","e","e","c","e","c"),
                 slope = round(c(beta_t, beta_tc, beta_g, beta_d, beta_dc, beta_p, beta_pc), 3),
                 r2 = round(c(r2_t, r2_tc, r2_g, r2_d, r2_dc, r2_p, r2_pc), 3))
res

# Plot
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(diag(gmat)), log10(diag(dmat)), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Dalechampia: Tulum CR", las=1, pch=16, col="blue3")
points(log10(var_g_g), log10(var_d_g))
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="firebrick")
legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors", "P eigenvectors"), 
       pch=c(16, 1, 16, 16), col=c("blue3", "black", "black", "firebrick"))

# Plot with modified axes
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(diag(gmat)), log10(diag(dmat)), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="", 
     ylab="",
     yaxt="n",
     xaxt="n",
     main="", las=1, pch=16, col="blue3")
mtext(expression(paste(italic(Dalechampia), " ", italic(scandens))), line=0.5, cex=.8)

points(log10(var_g_g), log10(var_d_g), pch=16)
points(log10(var_g_d), log10(var_d_d))
points(log10(var_g_p), log10(var_d_p), pch=16, col="firebrick")

mean1 = mean(log10(c(diag(dmat), var_d_g, var_d_d)))
mean2 = mean(log10(c(diag(gmat), var_g_g, var_g_d)))
segments(x0=mean2-10, y0=mean1-10, x1=mean2+10, y1=mean1+10)

legend("bottomright", legend=c(paste0("Original traits (", round(100*r2_t, 1),"%)"),
                               paste0("G directions (", round(100*r2_g, 1),"%)"),
                               paste0("D directions (", round(100*r2_d, 1),"%)"),
                               paste0("P directions (", round(100*r2_p, 1),"%)")),
       pch=c(16, 16, 1, 16), col=c("blue3", "black", "black", "firebrick"))

axis(1, at=c(-1, -.5, 0), signif(10^c(-1, -.5, 0),1))

#xt3 = c(1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
#x3at = log10(100*log(xt3)^2/(2/pi))
#axis(2, at=x3at, signif(xt3, 4), las=1)

x3at = seq(-1, 1, .5)
x3 = exp(sqrt(((10^x3at)/100)*(2/pi)))
axis(2, at=x3at, signif(x3, 3), las=1)

#exp(sqrt(diag(dmat/100)*(2/pi)))



# Cond
xmin = log10(min(c(var_g_g_c, var_g_d_c), na.rm=T))
xmax = log10(max(c(var_g_g_c, var_g_d_c), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(diag(gmat)), log10(diag(dmat)), col="white", pch=16, 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Conditional evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Dalechampia: Tulum CR", las=1)

legend("bottomright", c("Original traits", "G eigenvectors", "D eigenvectors", "P eigenvectors"), 
       pch=c(16, 1, 16, 16), col=c("blue3", "black", "black", "firebrick"))

points(log10(cvals), log10(diag(dmat)), pch=16, col="blue")
arrows(log10(diag(gmat)), log10(diag(dmat)), log10(cvals), log10(diag(dmat)), code=2, length=.05)


points(log10(var_g_d), log10(var_d_d), pch=16, col="grey")
points(log10(var_g_d_c), log10(var_d_d), pch=16)
arrows(log10(var_g_d), log10(var_d_d), log10(var_g_d_c), log10(var_d_d), code=2, length=.05)

points(log10(var_g_p), log10(var_d_p), pch=16, col="grey")
points(log10(var_g_p_c), log10(var_d_p), pch=16, col="firebrick")
arrows(log10(var_g_p), log10(var_d_p), log10(var_g_p_c), log10(var_d_p), code=2, length=.05)


points(log10(var_g_g_c), log10(var_d_g))




#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)


######## Analyse G vs. D: Mexico D. scandens ########

# The D matrix
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)

load(file="analyses/dalechampia/Dmat_Mexico_SG_75k.RData")
n=5
dpost = mod$VCV[, 1:(n*n)]
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("BA", "GA", "GSD", "SW", "GAD")
dmat

#plot(mod$VCV[,1])

#The G matrix
dat = read.csv("data/dalechampia/tovar/Merida diallel.csv", header=T)
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$Gwtot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tovar.RData")

n = 5
gpost = mod$VCV[, 1:(n*n)]
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD")
round(gmat, 2)

#The P matrix
load("analyses/dalechampia/MeanP_SG.RData")
pmat = MeanP

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)
evolvabilityMeans(pmat)
cov2cor(gmat)
cov2cor(dmat)
cov2cor(pmat)

vals = computeGD(gmat, dmat, pmat, nPop=12, nFam = 45, species="Dalechampia scandens B", SE=F, plot=F)

#Uncertainty over the posterior
out = list()
for(i in 1:100){
  sgmat = matrix(gpost[i,], nrow=n)
  sdmat = matrix(dpost[i,], nrow=n)
  #sdmat = dmat
  out[[i]] = computeGD(sgmat, sdmat, pmat, species="Dalechampia scandens A")   
}

slopes = lapply(out, function(x) x$res$slope)
slopemean = apply(simplify2array(slopes), 1, median)
slopeSE = apply(simplify2array(slopes), 1, sd)

vals$res$slope_MC = slopemean
vals$res$SE = slopeSE

vals

gdDF = data.frame(species="Dalechampia_scandens_B", g = "Dalechampia scandens B: Tovar", ntraits = ncol(gmat), 
                  dims = "line",
                  ndims = 1,
                  traitgroups = "flo",
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Dalechampia scandens B: Bolstad et al. Mexico greenhouse", nPop = 12, 
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

save(gdDF, file="analyses/dalechampia/gdDF_Tovar_MX.RData")

# Compute eigenvectors etc.
g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e

p_ev = eigen(pmat)$vectors
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
     ylab="log10 (Divergence [%])", 
     main="Dalechampia: Tovar Mexico", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

# Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)


#### Divergence vectors Tulum MX ####

# The G matrix
dat = read.csv("data/dalechampia/tulum/TulumDiallel.csv", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tulum.RData")
n = 6
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

names(dat)
dat=dat[,which(colnames(dat) %in% colnames(gmat))]
dat=dat[,c(6,5,3,4,1,2)]
names(dat)==colnames(gmat)
z0 = colMeans(dat, na.rm=T)

# Pop means
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$ind = paste(dat$pop, dat$plant, dat$seed, sep="_")
dat$ind = factor(dat$ind)
dat = dat[dat$population %in% c("BA", "C", "CA", "CC", "CN", "CP", "GS", "LM", "M", "PM", "T", "Tulum51"),]
dat$population = factor(dat$population)
head(dat)
means = apply(dat[,which(colnames(dat) %in% colnames(gmat))], 2, function(x) tapply(x, dat$population, mean, na.rm=T))
means = means[,c(6,5,3,4,1,2)]
head(means)
dim(means)

source("code/computeDelta.R")
outdat = computeDelta2(gmat/100, means, z0)

deltaDF = data.frame(species="Dalechampia_scandens_A", g="Dalechampia scandens A: Tulum", ntraits=ncol(gmat), 
                     d="Dalechampia scandens A: Bolstad et al. Mexico greenhouse", pop=rownames(means), 
                     emean=outdat$emean,
                     emin=outdat$emin,
                     emax=outdat$emax,
                     cmean=outdat$cmean,
                     div=outdat$div, edelta=outdat$edelta, cdelta=outdat$cdelta,
                     theta=outdat$theta, row.names=NULL)

head(deltaDF)
save(deltaDF, file="analyses/dalechampia/deltaDF_Tulum_MX.RData")

x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0,2), las=1,
     xlab="", ylab="",
     main="Dalechampia Tulum MX")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 56
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

#### Divergence vectors Tulum CR ####

# The G matrix
dat = read.csv("data/dalechampia/tulum/TulumDiallel.csv", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tulum.RData")
n = 6
gpost = mod$VCV[, 1:(n*n)]
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

names(dat)
dat = dat[,which(colnames(dat) %in% colnames(gmat))]
dat = dat[,c(6,5,3,4,1,2)]
names(dat)==colnames(gmat)
z0 = colMeans(dat, na.rm=T)

# Pop means
dat = read.table("data/dalechampia/populations/costarica_greenhouse.txt", header=T)
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GW)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$ind = dat$ID
dat$ind = factor(dat$ind)
dat = dat[dat$Pop %in% c("S1","S2","S3","S6","S7","S8","S9","S11","S12","S13","S16","S19","S20","S21","S22","S26"),]
dat$population = factor(dat$Pop)

means = apply(dat[,which(colnames(dat) %in% colnames(gmat))], 2, function(x) tapply(x, dat$population, mean, na.rm=T))
means = means[,c(6,5,3,4,2,1)]
colnames(means)==colnames(gmat)
head(means)
dim(means)

source("code/computeDelta.R")
outdat = computeDelta2(gmat/100, means, z0)

outlist = list()
for(i in 1:100){
  sgmat = matrix(gpost[i,], nrow=n)/100
  outlist[[i]] = computeDelta2(sgmat/100, means, z0)
}
meanmat = data.frame(apply(simplify2array(lapply(outlist, as.matrix)), 1:2, mean))
semat = data.frame(apply(simplify2array(lapply(outlist, as.matrix)), 1:2, sd))

rele = lapply(outlist, function(x) mean(x$edelta/x$emean))
mean(unlist(rele))
sd(unlist(rele))
relec = lapply(outlist, function(x) mean(x$cdelta/x$cmean))
mean(unlist(relec))
sd(unlist(relec))


deltaDF = data.frame(species="Dalechampia_scandens_A", g="Dalechampia scandens A: Tulum", ntraits=ncol(gmat), 
                     d="Dalechampia scandens A: Opedal et al. Costa Rica greenhouse", pop=rownames(means), 
                     emean=outdat$emean,
                     emin=outdat$emin,
                     emax=outdat$emax,
                     cmean=outdat$cmean,
                     div=outdat$div, edelta=outdat$edelta, cdelta=outdat$cdelta,
                     theta=outdat$theta, row.names=NULL)

head(deltaDF)
save(deltaDF, file="analyses/dalechampia/deltaDF_Tulum_CR.RData")

x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0,2), las=1,
     xlab="", ylab="",
     main="Dalechampia Tulum CR")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 38
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

#### Divergence vectors Tovar MX ####

# The G matrix
dat = read.csv("data/dalechampia/tovar/Merida diallel.csv", header=T)

# Define some traits
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$Gwtot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tovar.RData")
n = 5
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD")
round(gmat, 2)

names(dat)
dat = dat[,which(colnames(dat) %in% colnames(gmat))]
dat = dat[,c(5,4,2,3,1)]
names(dat)==colnames(gmat)
z0 = colMeans(dat, na.rm=T)

# Pop means
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

dat$ind = paste(dat$pop, dat$plant, dat$seed, sep="_")
dat$ind = factor(dat$ind)
dat = dat[dat$population %in% c("CH", "CO", "D", "E", "HO", "MO", "P", "S", "V", "X", "XU", "Tovar70"),]
dat$population = factor(dat$population)
head(dat)
means = apply(dat[,which(colnames(dat) %in% colnames(gmat))], 2, function(x) tapply(x, dat$population, mean, na.rm=T))
means = means[,c(5,4,2,3,1)]
head(means)
dim(means)

outdat = computeDelta2(gmat/100, means, z0)

deltaDF = data.frame(species="Dalechampia_scandens_B", g="Dalechampia scandens B: Tovar", ntraits=ncol(gmat), 
                     d="Dalechampia scandens B: Bolstad et al. Mexico greenhouse", pop=rownames(means), 
                     emean=evolvabilityMeans(gmat)[1],
                     emin=evolvabilityMeans(gmat)[2],
                     emax=evolvabilityMeans(gmat)[3],
                     cmean=evolvabilityMeans(gmat)[4],
                     div=outdat[,1], edelta=outdat[,2], cdelta=outdat[,3], 
                     theta=outdat[,4], row.names=NULL)
head(deltaDF)

save(deltaDF, file="analyses/dalechampia/deltaDF_Tovar_MX.RData")

x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0,2), las=1,
     xlab="", ylab="",
     main="Dalechampia Tovar MX")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 24
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")


#### Selection/Evolvability analysis ####
# The D matrix
dat = read.table("data/dalechampia/populations/costarica_greenhouse.txt", header=T)
dat$GSD = apply(subset(dat, select=c(GSDl, GSDc, GSDr)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHl, GHr)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWl, SWc, SWr)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GW)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)

load(file="analyses/dalechampia/Dmat_CostaRica_75k.RData")
n=6
dpost = mod$VCV[,1:(n*n)]
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
dmat

#plot(mod$VCV[,1])

# The G matrix
dat = read.csv("data/dalechampia/tulum/TulumDiallel.csv", header=T)
dat$GSD = apply(subset(dat, select=c(GSDL, GSDC, GSDR)), 1, mean, na.rm=T)
dat$GH  = apply(subset(dat, select=c(GHL, GHR)), 1, mean, na.rm=T)
dat$SW  = apply(subset(dat, select=c(SWL, SWC, SWR)), 1, mean, na.rm=T)
dat$UBL = apply(subset(dat, select=c(UBL, UBC, UBR)), 1, mean, na.rm=T)
dat$LBL = apply(subset(dat, select=c(LBL, LBC, LBR)), 1, mean, na.rm=T)
dat$GA  = sqrt(dat$GH*dat$GWTot)
dat$BA  = with(dat, sqrt(UBL*UBW)+sqrt(LBL*LBW))

load("analyses/dalechampia/Gmat_Tulum.RData")

n = 6
gpost = mod$VCV[,1:(n*n)]
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

# The mean P matrix
load(file="analyses/dalechampia/MeanP.RData")
pmat = MeanP

# Selection data
flosel = read.table("C:/data/FlowerSelection/flowerbeta.txt", header=T)
#flosel = flosel[flosel$Study!="Soteras_et_al_2020",]
flosel$StudySpecies = as.factor(paste(flosel$Study, flosel$Genus, flosel$epithet, sep="_"))
flosel$BetaRaw = flosel$BetaVar/flosel$SD
flosel$BetaRawSE = flosel$BetaVarSE/flosel$SD
flosel$CV = flosel$SD/flosel$Mean
ww = which(flosel$Study=="Albertsen_et_al._201x")
flosel$PL[ww] = 1-((0.13*flosel$pollen[ww])/(1+(0.13*flosel$pollen[ww])))
sa = flosel[flosel$Study=="Albertsen_et_al._201x",]

sa$fit = sa$PollSize-sa$Mean

names(sa)

bdat = subset(sa, select=c("Population", "Year", "Trait", "BetaMean"))
head(bdat, 9)

bdat = data.frame(matrix(bdat$BetaMean, nrow=8))
colnames(bdat)=unique(sa$Trait)
bdat

Bmat = cov(bdat)
round(Bmat, 3)

names(bdat)
Gred = gmat[c(1:3, 6), c(1:3, 6)]


Rmat = Gred %*% Bmat %*% Gred
Rmat

# Derive eigenvectors
# Compute eigenvectors etc.
gmat = Gred
dmat = dmat[c(1:3, 6), c(1:3, 6)]
bmat = Bmat
rmat = Rmat
pmat = pmat[c(1:3, 6), c(1:3, 6)]

library(evolvability)

g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
var_r_g = evolvabilityBeta(rmat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e
var_r_d = evolvabilityBeta(rmat, Beta = d_ev)$e

r_ev = eigen(rmat)$vectors
var_g_r = evolvabilityBeta(gmat, Beta = r_ev)$e
var_d_r = evolvabilityBeta(dmat, Beta = r_ev)$e
var_r_r = evolvabilityBeta(rmat, Beta = r_ev)$e

p_ev = eigen(pmat)$vectors
var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e
var_r_p = evolvabilityBeta(rmat, Beta = p_ev)$e

# Plot R vs. D

plot(log10(var_r_r), log10(var_d_r), xlim=c(-4, 0), ylim=c(-1, 2))
points(log10(var_r_d), log10(var_d_d), pch=16)
points(log10(var_r_p), log10(var_d_p), pch=16, col="firebrick")
points(log10(diag(rmat)), log10(diag(dmat)), pch=16, col="blue")

x11(width=5, height=5)
xmin = log10(min(c(var_r_r, var_r_d), na.rm=T))
xmax = log10(max(c(var_r_r, var_r_d), na.rm=T))
ymin = log10(min(c(var_d_r, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_r, var_d_d), na.rm=T))
plot(log10(var_r_r), log10(var_d_r), 
     xlim=c(xmin, xmax), ylim=c(ymin-1, ymax), 
     xlab="log10 (Exp. divergence)", 
     ylab="log10 (Divergence [x100])", 
     main="Dalechampia: Tulum CR", las=1)
points(log10(var_r_d), log10(var_d_d), pch=16)
points(log10(var_r_p), log10(var_d_p), pch=16, col="firebrick")
points(log10(diag(rmat)), log10(diag(dmat)), pch=16, col="blue3")

legend("bottomright", c("Original traits", "R eigenvectors", "D eigenvectors", "P eigenvectors"), 
       pch=c(16, 1, 16, 16), col=c("blue3", "black", "black", "firebrick"))



var_b_p = evolvabilityBeta(bmat, Beta = p_ev)$e

plot(log10(var_g_p), log10(var_b_p), xlim=c(-1.5,.5), ylim=c(-2,0))
points(log10(diag(gmat)), log10(diag(bmat)), pch=16)

x11(width=12*.75, height=4.5*.75)
par(mfrow=c(1,3))
plot(log10(var_g_p), log10(var_d_p), las=1, ylim=c(-.3, 1), xlab="", ylab="")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue3")
mtext("log10 (Evolvability)", 1, line=2.5)
mtext("log10 (Divergence)", 2, line=2.5)

plot(log10(var_b_p), log10(var_d_p), las=1, ylim=c(-.3, 1), xlab="", ylab="")
points(log10(diag(bmat)), log10(diag(dmat)), pch=16, col="blue3")
mtext("log10 (Var in selection)", 1, line=2.5)

plot(log10(var_r_p), log10(var_d_p), las=1, ylim=c(-.3, 1), xlab="", ylab="")
points(log10(diag(rmat)), log10(diag(dmat)), pch=16, col="blue3")
mtext("log10 (Pred. divergence)", 1, line=2.5)


mt = lm(log(diag(dmat))~log(diag(gmat)))
beta_t = summary(mt)$coef[2,1]
beta_t
r2_t = summary(mt)$r.squared
r2_t

mtc = lm(log(diag(dmat))~log(cvals))
beta_tc = summary(mtc)$coef[2,1]
beta_tc
r2_tc = summary(mtc)$r.squared
r2_tc

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

mdc = lm(log(var_d_d)~log(var_g_d_c))
beta_dc = summary(mdc)$coef[2,1]
beta_dc
r2_dc = summary(mdc)$r.squared
r2_dc

mp = lm(log(var_d_p)~log(var_g_p))
beta_p = summary(mp)$coef[2,1]
beta_p
r2_p = summary(mp)$r.squared
r2_p

mpc = lm(log(var_d_p)~log(var_g_p_c))
beta_pc = summary(mpc)$coef[2,1]
beta_pc
r2_pc = summary(mpc)$r.squared
r2_pc

res = data.frame(traits = c("Original", "Original","G eigenvectors", "D eigenvectors", "D eigenvectors", "P eigenvectors", "P eigenvectors"),
                 evol = c("e","c","e","e","c","e","c"),
                 slope = round(c(beta_t, beta_tc, beta_g, beta_d, beta_dc, beta_p, beta_pc), 3),
                 r2 = round(c(r2_t, r2_tc, r2_g, r2_d, r2_dc, r2_p, r2_pc), 3))
res
