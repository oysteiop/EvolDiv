##############################################
#### Dalechampia 'scandensoides' analyses ####
##############################################

rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)

########################################################
###### G matrix for the Tulum population (Mexico) ######
########################################################
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

###########################################################
###### G matrix for the Tovar population (Venezuela) ######
###########################################################
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



##################################################################
###### D matrix for Mexican populations of D. scandensoides ######
##################################################################
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

#############################################################
###### D matrix for Mexican populations of D. scandens ######
#############################################################
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

########################################################
###### D matrix for Costa Rican populations ######
########################################################
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

#########################################
######## Analyse G vs. D: Mexico ########
#########################################

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

n = length(nms)
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)/100
colnames(gmat) = rownames(gmat) = c("BA", "GA", "GSD", "SW", "GAD", "ASD")
round(gmat, 2)

pmatdat = subset(dat, select=c("BA", "GA", "GSD", "SW", "GAD", "ASD"))
#pmatdat = apply(pmatdat, 2, function(x) log(x))
pmat = cov(pmatdat)

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)

#Compute eigenvectors etc.
g_ev=eigen(gmat)$vectors
var_g_g=evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
#var_d_g = diag(t(g_ev) %*% dmat %*% g_ev)

d_ev=eigen(dmat)$vectors
var_g_d=evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d=evolvabilityBeta(dmat, Beta = d_ev)$e
#var_d_d = diag(t(d_ev) %*% dmat %*% d_ev)

p_ev=eigen(pmat)$vectors
var_g_p=evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p=evolvabilityBeta(dmat, Beta = p_ev)$e

#Compute summary stats
mg=lm(log(var_d_g)~log(var_g_g))
beta_g=summary(mg)$coef[2,1]
beta_g
r2_g=summary(mg)$r.squared
r2_g

md=lm(log(var_d_d)~log(var_g_d))
beta_d=summary(md)$coef[2,1]
beta_d
r2_d=summary(md)$r.squared
r2_d

mp=lm(log(var_d_p)~log(var_g_p))
beta_p=summary(mp)$coef[2,1]
beta_p
r2_p=summary(mp)$r.squared
r2_p

#Plot both relationship,s upper/lower bounds of scaling relationship
#Plot
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
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

points(log10(var_g_p), log10(var_d_p), pch=16, col="red")

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

#############################################
######## Analyse G vs. D: Costa Rica ########
#############################################

#The D matrix
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

#Compute eigenvectors etc.
g_ev=eigen(gmat)$vectors
var_g_g=evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
#var_d_g = diag(t(g_ev) %*% dmat %*% g_ev)

d_ev=eigen(dmat)$vectors
var_g_d=evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d=evolvabilityBeta(dmat, Beta = d_ev)$e
#var_d_d = diag(t(d_ev) %*% dmat %*% d_ev)

p_ev=eigen(pmat)$vectors
var_g_p=evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p=evolvabilityBeta(dmat, Beta = p_ev)$e

#Compute summary stats
mg=lm(log(var_d_g)~log(var_g_g))
beta_g=summary(mg)$coef[2,1]
beta_g
r2_g=summary(mg)$r.squared
r2_g

md=lm(log(var_d_d)~log(var_g_d))
beta_d=summary(md)$coef[2,1]
beta_d
r2_d=summary(md)$r.squared
r2_d

mp=lm(log(var_d_p)~log(var_g_p))
beta_p=summary(mp)$coef[2,1]
beta_p
r2_p=summary(mp)$r.squared
r2_p

#Plot both relationship,s upper/lower bounds of scaling relationship
#Plot
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
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

points(log10(var_g_p), log10(var_d_p), pch=16, col="red")

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

acos(t(g_ev[,2]) %*% d_ev[,2])*(180/pi)

# Ellipse plot
source(GDellipse)
GDellipse(dmat, gmat, xlim=c(-2,2), ylim=c(-1.5, 1.5), main="Dalechampia: Tulum: Costa Rica")



#####################################################
######## Analyse G vs. D: Mexico D. scandens ########
#####################################################

#The D matrix
dat = read.table("data/dalechampia/populations/mexico_greenhouse.txt", header=T)

load(file="analyses/dalechampia/Dmat_Mexico_SG_75k.RData")
n=5
dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = c("BA", "GA", "GSD", "SW", "GAD")
dmat

plot(mod$VCV[,1])

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

#Compute eigenvectors etc.
g_ev=eigen(gmat)$vectors
var_g_g=evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
#var_d_g = diag(t(g_ev) %*% dmat %*% g_ev)

d_ev=eigen(dmat)$vectors
var_g_d=evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d=evolvabilityBeta(dmat, Beta = d_ev)$e
#var_d_d = diag(t(d_ev) %*% dmat %*% d_ev)

p_ev=eigen(pmat)$vectors
var_g_p=evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p=evolvabilityBeta(dmat, Beta = p_ev)$e

#Compute summary stats
mg=lm(log(var_d_g)~log(var_g_g))
beta_g=summary(mg)$coef[2,1]
beta_g
r2_g=summary(mg)$r.squared
r2_g

md=lm(log(var_d_d)~log(var_g_d))
beta_d=summary(md)$coef[2,1]
beta_d
r2_d=summary(md)$r.squared
r2_d

mp=lm(log(var_d_p)~log(var_g_p))
beta_p=summary(mp)$coef[2,1]
beta_p
r2_p=summary(mp)$r.squared
r2_p

#Plot both relationship,s upper/lower bounds of scaling relationship
#Plot
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
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

points(log10(var_g_p), log10(var_d_p), pch=16, col="red")

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
