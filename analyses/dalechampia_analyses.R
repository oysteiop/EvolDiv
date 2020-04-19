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

###### D matrix for Mexican populations of D. scandensoides ######
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

# D matrix
unique(dat$population)
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

# D matrix
unique(dat$population)
dat = dat[dat$population %in% c("CH", "CO", "D", "E", "HO", "MO", "P", "S", "V", "X", "XU", "Tovar70"),]

X = subset(dat, select=c(BA, GA, GSD, SW, GAD, population, ind))
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
mod<-MCMCglmm(cbind(BA, GA, GSD, SW, GAD) ~ -1+trait,
              random= ~us(trait):population + us(trait):ind,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/dalechampia/Dmat_Mexico_SG_75k.RData")

###### D matrix for Costa Rican populations ######
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

# D matrix
unique(dat$population)
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
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
dmat

plot(mod$VCV[,1])

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
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

pmatdat = subset(dat, select=c("BA", "GA", "GSD", "SW", "GAD", "ASD"))
#pmatdat = apply(pmatdat, 2, function(x) log(x))
pmat = cov(pmatdat)

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)

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

source("code/plot_GD.R")
vals = plot_GD(gmat, dmat, species="Dalechampia scandens A", plot=T)

gdDF = data.frame(sp="Dalechampia_scandens_A", g = "Dalechampia scandens: Tulum", traits = ncol(gmat), 
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Dalechampia scandens A: MX", npops = 12, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaG = vals[1,1], r2G = vals[2,1],
                  betaD = vals[3,1], r2D = vals[4,1],
                  theta = vals[5,1], row.names = NULL)
head(gdDF)

save(gdDF, file="analyses/dalechampia/gdDF_Tulum_MX.RData")

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
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

pmatdat = subset(dat, select=c("BA", "GA", "GSD", "SW", "GAD", "ASD"))
#pmatdat = apply(pmatdat, 2, function(x) log(x))
pmat = cov(pmatdat)

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)

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
     main="Dalechampia: Tulum CR", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

cvals=NULL
for(i in 1:ncol(gmat)){
        b=rep(0,ncol(gmat))
        b[i]=1
        cvals[i] = evolvabilityBeta(gmat, b)$c
}
points(log10(cvals), log10(diag(dmat)), pch=16, col="red")

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

vals = plot_GD(gmat, dmat, species="Dalechampia scandens A", plot=T)

gdDF = data.frame(sp="Dalechampia_scandens_A", g = "Dalechampia scandens: Tulum", traits = ncol(gmat), 
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Dalechampia scandens A: CR", npops = 17, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaG = vals[1,1], r2G = vals[2,1],
                  betaD = vals[3,1], r2D = vals[4,1],
                  theta = vals[5,1], row.names = NULL)
head(gdDF)

save(gdDF, file="analyses/dalechampia/gdDF_Tulum_CR.RData")

######## Analyse G vs. D: Mexico D. scandens ########

# The D matrix
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)

load(file="analyses/dalechampia/Dmat_Mexico_SG_75k.RData")
n=5
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
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD")
round(gmat, 2)

pmatdat = na.omit(subset(dat, select=c("BA", "GA", "GSD", "SW", "GAD")))
#pmatdat = apply(pmatdat, 2, function(x) log(x))
pmat = cov(pmatdat)

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)
cov2cor(gmat)
cov2cor(dmat)

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

vals = plot_GD(gmat, dmat, species="Dalechampia scandens B", plot=T)

gdDF = data.frame(sp="Dalechampia_scandens_B", g = "Dalechampia scandens: Tovar", traits = ncol(gmat), 
                  emean = evolvabilityMeans(gmat)[1],
                  emin = evolvabilityMeans(gmat)[2],
                  emax = evolvabilityMeans(gmat)[3],
                  cmean = evolvabilityMeans(gmat)[4],
                  imean = evolvabilityMeans(gmat)[7],
                  d = "Dalechampia scandens B: MX", npops = 12, 
                  dmean = evolvabilityMeans(dmat)[1],
                  betaG = vals[1,1], r2G = vals[2,1],
                  betaD = vals[3,1], r2D = vals[4,1],
                  theta = vals[5,1], row.names = NULL)
head(gdDF)

save(gdDF, file="analyses/dalechampia/gdDF_Tovar_MX.RData")

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

outdat = computeDelta(gmat/100, means, z0)

deltaDF = data.frame(sp="Dalechampia_scandens_A", g="Dalechampia scandens: Tulum", traits=ncol(gmat), 
                     d="Dalechampia scandens A: MX", pop=rownames(means), 
                     emean=evolvabilityMeans(gmat)[1],
                     emin=evolvabilityMeans(gmat)[2],
                     emax=evolvabilityMeans(gmat)[3],
                     cmean=evolvabilityMeans(gmat)[4],
                     div=outdat[,1], edelta=outdat[,2], cdelta=outdat[,3], 
                     theta=outdat[,4], row.names=NULL)
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

outdat = computeDelta(gmat/100, means, z0)

deltaDF = data.frame(sp="Dalechampia_scandens_A", g="Dalechampia scandens: Tulum", traits=ncol(gmat), 
                     d="Dalechampia scandens A: CR", pop=rownames(means), 
                     emean=evolvabilityMeans(gmat)[1],
                     emin=evolvabilityMeans(gmat)[2],
                     emax=evolvabilityMeans(gmat)[3],
                     cmean=evolvabilityMeans(gmat)[4],
                     div=outdat[,1], edelta=outdat[,2], cdelta=outdat[,3], 
                     theta=outdat[,4], row.names=NULL)
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

outdat = computeDelta(gmat/100, means, z0)

deltaDF = data.frame(sp="Dalechampia_scandens_B", g="Dalechampia scandens: Tovar", traits=ncol(gmat), 
                     d="Dalechampia scandens B: MX", pop=rownames(means), 
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
