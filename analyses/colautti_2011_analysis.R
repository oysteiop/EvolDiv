##############################################################
#### Lythrum salicaria: Data from Colautti & Barrett 2011 ####
##############################################################

rm(list=ls())
library(evolvability)
library(plyr)
library(reshape2)
library(MCMCglmm)

# Summary stats
#Gdat = read.table("data/colautti/AllDataFixed.txt", dec=".", sep=",", header=T)
Gdat = read.table("data/colautti/AllDataFixed2.txt", header=T)
names(Gdat)

popmeans = as.data.frame(apply(Gdat[,c(11:23, 25:38)], 2, function(x) tapply(x, Gdat$Pop, mean, na.rm = T)))
head(popmeans)

popmeans = subset(popmeans, select=c("TLeafArea", "THeight", "Height2wk", "Height4wk",
                                     "FDays", "FStemWidth", "FVeg",
                                     "FInf", "HVeg","HInf","HVegW","HInfW"))
head(t(popmeans))
#write.csv2(melt(t(popmeans)), file = "data/colautti/popmeans.csv", row.names=T)

popvars = as.data.frame(apply(Gdat[,c(11:23, 25:38)], 2, function(x) tapply(x, Gdat$Pop, var, na.rm = T)))
head(popvars)

popvars = subset(popvars, select=c("TLeafArea", "THeight", "Height2wk", "Height4wk",
                                   "FDays", "FStemWidth", "FVeg",
                                   "FInf", "HVeg","HInf","HVegW","HInfW"))
head(t(popvars))
#write.csv2(melt(t(popvars)), file = "data/colautti/popvars.csv", row.names=T)

popn = as.data.frame(apply(Gdat[,c(11:23, 25:38)], 2, function(x) tapply(x>0, Gdat$Pop, sum, na.rm = T)))
head(popn)

popn = subset(popn, select=c("TLeafArea", "THeight", "Height2wk", "Height4wk",
                             "FDays", "FStemWidth", "FVeg",
                             "FInf", "HVeg","HInf","HVegW","HInfW"))
head(t(popn))
#write.csv2(melt(t(popn)), file = "data/colautti/popn.csv", row.names=T)

means = apply(Gdat[,c(11:23, 25:38)], 2, mean, na.rm=T)
means
cbind(means[c(23,1,5,7,9,14,11,12,16,17,19,20)])

vars = apply(Gdat[,c(11:23, 25:38)], 2, var, na.rm=T)
round(cbind(vars[c(23,1,5,7,9,14,11,12,16,17,19,20,26,27)]),2)

ev = c((vars[23]*0.088*4)/(means[23]^2)*100, #Leaf area at transplant, area
       (vars[1]*0.109*4)/(means[1]^2)*100, #Transplant height, linear
       (vars[5]*0.069*4)/(means[5]^2)*100, #Height at 2 weeks, linear
       (vars[7]*0.068*4)/(means[7]^2)*100, #Height at 4 weeks, linear
       (vars[9]*0.129*4)/(means[9]^2)*100, #Days to flowering, count
       (vars[14]*0.100*4)/(means[14]^2)*100, #Stem width at maturity, linear
       (vars[11]*0.111*4)/(means[11]^2)*100, #Veg size at maturity, linear
       (vars[12]*0.044*4)/(means[12]^2)*100, #Inflorescence size at maturity,linear
       (vars[16]*0.104*4)/(means[16]^2)*100, #Veg size at harvest, linear
       (vars[17]*0.060*4)/(means[17]^2)*100, #Inflorescence size at harvest, linear
       (vars[26]*0.070*4)*100, #Final vegetative biomass, mass_volume
       (vars[27]*0.057*4)*100) #Final inflorescence biomass, mass_volume

ev

#### Genetic variance for one population #### 
library(lme4)
one = subset(Gdat, Pop=="ONTI")
one = droplevels(one)
one$cID = 1:nrow(one)

names(one)
mod = lmer(THeight ~ Table + (1|TrueFam), data=one)
summary(mod)
Vg = VarCorr(mod)$TrueFam[1]*2
mu = mean(one$THeight, na.rm=T)
Vg/(mu^2)*100

# All populations
pops = unique(Gdat$Pop)
traits = c("TLeafLength", "TLeafWidth","THeight", "Height2wk", "Height4wk",
           "FDays", "FStemWidth", "FVeg", "FInf", "HVeg","HInf")

emat = matrix(NA, nrow=20, ncol=11)

for(p in 1:20){
  pop=pops[p]
    for(t in 1:11){
      tr=traits[t]  
      one = subset(Gdat, Pop==pop)
      one = droplevels(one)
      
      names(one)
      w=which(colnames(one)==tr)
      mod = lmer(one[,w] ~ Table + (1|TrueFam), data=one)
      Vg=VarCorr(mod)$TrueFam[1]*2
      mu=mean(one$THeight, na.rm=T)
      emat[p,t]=Vg/(mu^2)*100
      }
}

rownames(emat) = pops
colnames(emat) = traits

evals = colMeans(emat)
evals

########################################################
###### Joint estimation of D and average G ######
########################################################
Gdat = read.table("data/colautti/AllDataFixed2.txt", header=T)
Gdat$TLeafArea = sqrt(Gdat$TLeafArea)
names(Gdat)
traits = c("TLeafArea","THeight", "Height2wk", "Height4wk",
           "FStemWidth", "FVeg", "FInf")

X = subset(Gdat, select=c(traits, "Pop", "TrueFam"))
X[,1:length(traits)] = apply(X[,1:length(traits)], 2, function(x) 10*log(x))
X = na.omit(X)
head(X)
str(X)

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(traits)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(TLeafArea,THeight,Height2wk,Height4wk,FStemWidth,FVeg,FInf) ~ -1+trait,
              random= ~us(trait):Pop + us(trait):TrueFam,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/colautti/mod75.RData")

#The D matrix
Gdat = read.table("data/colautti/AllDataFixed2.txt", header=T)
names(Gdat)

traits = c("TLeafArea","THeight", "Height2wk", "Height4wk",
           "FStemWidth", "FVeg", "FInf")
n = length(traits)

load(file="analyses/colautti/mod75.RData")

dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = traits
dmat

# The (within-pop mean) G matrix
gmat = matrix(apply(mod$VCV, 2, median)[((n*n)+1):((n*n)*2)], nrow=n)
colnames(gmat) = rownames(gmat) = traits
round(gmat, 2)

# The (within-individual) residual matrix
rmat = matrix(apply(mod$VCV, 2, median)[(((n*n)*2)+1):((n*n)*3)], nrow=n)
colnames(rmat) = rownames(rmat) = traits
round(rmat, 2)

diag(gmat)/(diag(dmat) + diag(gmat))

pmatdat = na.omit(subset(Gdat, select=c(traits)))
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
     main="Lythrum salicaria", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

points(log10(var_g_p), log10(var_d_p), pch=16, col="red")

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

# Ellipse plot
source("GDellipse.R")
GDellipse(dmat, gmat, xlim=c(-2,2), ylim=c(-2, 2), main="Lythrum salicaria")


