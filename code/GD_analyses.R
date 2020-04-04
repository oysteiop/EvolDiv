#######################################
#### Analyses of G and D matrices #####
#######################################

rm(list=ls())
library(reshape2)
library(MCMCglmm)
library(evolvability)
source("code/computeGD.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")

# Species present in both databases
gsp = unlist(lapply(EVOBASE, function(x) x$Species))
dsp = unlist(lapply(POPBASE, function(x) x$Species))
both_sp = unique(gsp[which(gsp %in% dsp)])
both_sp

POPBASE = POPBASE[-c(14:15)] #Dropping second M. guttatus study

out = computeGD(species = both_sp[1], gmatrix = 1, dmatrix = 3)
out

#### Lobelia ####
out = computeGD(species = both_sp[1], gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Mean G-matrix
glist = list()
glist[[1]] = computeGD(species = both_sp[1], gmatrix = 1, dmatrix = 2)$G
glist[[2]] = computeGD(species = both_sp[1], gmatrix = 2, dmatrix = 2)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
unlist(lapply(POPBASE, function(x)x$Study_ID))
s=25 #Lobelia, Caruso 2012 Greenhouse
s=21 #Caruso 2003 Field 

means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means), colnames(out$D)==names(eV))

source("code/estimateD.R")
modD = estimateD(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Lobelia.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = gmat*100
out$G = out$G*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G),2)
signif(cov2cor(out$D),2)

# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

# Plot
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Lobelia siphilitica: Caruso 2003 D", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#### Brassica ####
out = computeGD(species = both_sp[13], gmatrix = 5, dmatrix = 1)

# Mean G matrix
glist = list()
for(i in 1:5){
  glist[[i]] = computeGD(species = both_sp[13], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
unlist(lapply(POPBASE, function(x)x$Study_ID))
s=30 #Brassica

means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

source("code/estimateD.R")
modD = estimateD(means, eV, thin=1000)

save(modD, file="analyses/adj_Dmats/Brassica.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Brassica.RData")

round(modD, 3)
eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = gmat*100
out$G = out$G*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
cov2cor(out$G)
cov2cor(out$D)

# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

# Plot
x11(width=5, height=5)
xmin=-1
#xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin=-2
#ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Brassica cretica", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#### Spergularia ####
out = computeGD(species = both_sp[7], gmatrix = 4, dmatrix = 1)
out$gmat

# Mean G matrix
gsp
EVOBASE[[19]]

glist = list()
for(i in 1:4){
  glist[[i]] = computeGD(species = both_sp[7], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
unlist(lapply(POPBASE, function(x)x$Study_ID))
s=8

means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

source("code/estimateD.R")
modD = estimateD(means, eV, thin=500)

save(modD, file="analyses/adj_Dmats/Spergularia.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Spergularia.RData")

round(modD, 3)

eigen(modD*100)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = gmat*100

out$G = out$G*100
#out$D = out$D*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
cov2cor(out$G)
cov2cor(out$D)

# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

# Plot
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin-2, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Spergularia marina", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#### Solanum ####
out = computeGD(species = both_sp[5], gmatrix = 3, dmatrix = 1)

# Mean G matrix
glist = list()
for(i in 1:3){
  glist[[i]] = computeGD(species = both_sp[5], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
unlist(lapply(POPBASE, function(x)x$Study_ID))
s=5

means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

source("code/estimateD.R")
modD = estimateD(means, eV, thin=1000)

save(modD, file="analyses/adj_Dmats/Solanum.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Solanum.RData")

round(modD, 3)

eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = gmat*100
out$G = out$G*100
#out$D = out$D*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
cov2cor(out$G)
cov2cor(out$D)

# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

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
     main="Solanum carolinense", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#### Clarkia ####
both_sp
out = computeGD(species = both_sp[9], gmatrix = 2, dmatrix = 1)
out$gmat

# Mean G matrix
glist = list()
for(i in 1:2){
  glist[[i]] = computeGD(species = both_sp[9], gmatrix = i, dmatrix = 1)$G
}
gmat = apply(simplify2array(glist), 1:2, mean)

# Estimate error-corrected D matrix
unlist(lapply(POPBASE, function(x)x$Study_ID))
s=9

means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

source("code/estimateD.R")
modD = estimateD(means, eV, thin=1000)

save(modD, file="analyses/adj_Dmats/Clarkia.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Clarkia.RData")

round(modD, 3)

eigen(modD)$values
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = gmat*100

out$G = out$G*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

# Plot
x11(width=5, height=5)
xmin=-1
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin=-2
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Clarkia dudleyana", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#### Ipomopsis ####
both_sp
out = computeGD(species = both_sp[6], gmatrix = 1, dmatrix = 4)

# Estimate error-corrected D matrix
unlist(lapply(POPBASE, function(x)x$Study_ID))
s=22 #Caruso 2000 field
s=23 #Caruso 2001 field
s=36 #Campbell unpublised 2019 field

means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
eV = POPBASE[[s]]$eV[,-1]
eV = eV[,match(colnames(out$D), names(eV))]
c(colnames(out$D)==names(means),colnames(out$D)==names(eV))

source("code/estimateD.R")
modD = estimateD(means, eV, thin=1000)

save(modD, file="analyses/adj_Dmats/Ipomopsis_Campbell2019.RData")

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2000.RData")
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2001.RData")
load(file="analyses/adj_Dmats/Ipomopsis_Caruso2019.RData")

round(modD, 3)

eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

out$D = modD*100
out$G = out$G*100
#out$D = out$D*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)
signif(cov2cor(out$G), 2)
signif(cov2cor(out$D), 2)

# Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g), na=na.exclude)
beta_g = summary(mg)$coef[2,1]
beta_g
r2_g = summary(mg)$r.squared
r2_g

md = lm(log(var_d_d)~log(var_g_d), na=na.exclude)
beta_d = summary(md)$coef[2,1]
beta_d
r2_d = summary(md)$r.squared
r2_d

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
     main="Ipomopsis aggregata: Campbell 2019 D", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

#### Evolvability along the vector of divergence ####

##### Mimulus guttatus ####
gsp

m1 = abs(EVOBASE[[3]]$Means[c(2,4:9)])
m2 = abs(EVOBASE[[4]]$Means[c(2,4:9)])
names(m1)==names(m2)

g1 = droptraits(EVOBASE[[3]]$G)
g2 = droptraits(EVOBASE[[4]]$G)
g1 = meanStdG(g1, m1)*100
g2 = meanStdG(g2, m2)*100
c(names(m1)==colnames(g1), names(m2)==colnames(g2))

delta = m1-m2
delta = delta/sqrt(sum(delta^2)) #Unit-length

evolvabilityBeta(g1, delta)$e
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[1]

evolvabilityBeta(g1, delta)$c
evolvabilityBeta(g1, delta)$c/evolvabilityMeans(g1)[4]

evolvabilityBeta(g2, delta)$e
evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[1]

evolvabilityBeta(g2, delta)$c
evolvabilityBeta(g2, delta)$c/evolvabilityMeans(g2)[4]

mG = apply(simplify2array(list(g1, g2)), 1:2, mean)
evolvabilityBeta(mG, delta)$e
evolvabilityMeans(mG)
evolvabilityBeta(mG, delta)$e/evolvabilityMeans(mG)[1]

plot(1:2, c(evolvabilityBeta(g1, delta)$e, evolvabilityBeta(g2, delta)$e), 
     pch=16, ylim=c(-3,30), xlim=c(0,10))
points(1:2, c(evolvabilityMeans(g1)[1], evolvabilityMeans(g2)[1]))
points(1:2, c(evolvabilityMeans(g1)[2], evolvabilityMeans(g2)[2]))
points(1:2, c(evolvabilityMeans(g1)[3], evolvabilityMeans(g2)[3]))

180-acos(t(eigen(g1)$vectors[,1]) %*% eigen(g2)$vectors[,1])*(180/pi) #Angle between the 2 G matrices
acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector
180-acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G2 vs. divergence vector

##### Mimulus micranthus ####
gsp
m1 = EVOBASE[[5]]$Means[c(2,4:9)]
m2 = EVOBASE[[6]]$Means[c(2,4:9)]
names(m1)==names(m2)

g1 = droptraits(EVOBASE[[3]]$G)
g2 = droptraits(EVOBASE[[4]]$G)
g1 = meanStdG(g1, m1)*100
g2 = meanStdG(g2, m2)*100

c(names(m1)==colnames(g1), names(m2)==colnames(g2))

delta = m2-m1
delta = delta/sqrt(sum(delta^2)) #Unit-length

evolvabilityBeta(g1, delta)$e
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[1]

evolvabilityBeta(g1, delta)$c

evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[1]

evolvabilityBeta(g2, delta)$c

points(3:4, c(evolvabilityBeta(g1, delta)$e, evolvabilityBeta(g2, delta)$e), pch=16)
points(3:4, c(evolvabilityMeans(g1)[1], evolvabilityMeans(g2)[1]))
points(3:4, c(evolvabilityMeans(g1)[2], evolvabilityMeans(g2)[2]))
points(3:4, c(evolvabilityMeans(g1)[3], evolvabilityMeans(g2)[3]))

180-acos(t(eigen(g1)$vectors[,1]) %*% eigen(g2)$vectors[,1])*(180/pi) #Angle between the 2 G matrices
acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector
180-acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G2 vs. divergence vector

##### Holcus lanatus ####
gsp
m1 = EVOBASE[[32]]$Means[c(2,5,6,7)]
m2 = EVOBASE[[33]]$Means[c(1,4,5,6)]
names(m1)==names(m2)

g1 = droptraits(EVOBASE[[32]]$G)
g1 = meanStdG(g1, m1)*100
c(names(m1)==colnames(g1))

delta = m2-m1
delta = delta/sqrt(sum(delta^2)) #Unit-length

evolvabilityBeta(g1, delta)$e
evolvabilityBeta(g1, delta)$c
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[1]

points(5, evolvabilityBeta(g1, delta)$e, pch=16)
points(5, evolvabilityMeans(g1)[1])
points(5, evolvabilityMeans(g1)[2])
points(5, evolvabilityMeans(g1)[3])

acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector

gsp
m1 = EVOBASE[[32]]$Means[c(3,5,6,7)]
m2 = EVOBASE[[33]]$Means[c(2,4,5,6)]
names(m1)==names(m2)

g2 = droptraits(EVOBASE[[33]]$G)
g2 = meanStdG(g2, m2)*100
c(names(m1)==colnames(g2))

delta = m2-m1
delta = delta/sqrt(sum(delta^2)) #Unit-length

evolvabilityBeta(g2, delta)$e
evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[1]

points(6, evolvabilityBeta(g2, delta)$e, pch=16)
points(6, evolvabilityMeans(g2)[1])
points(6, evolvabilityMeans(g2)[2])
points(6, evolvabilityMeans(g2)[3])

acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector

##### Nigella degenii ####
g1 = droptraits(EVOBASE[[23]]$G)

m1 = EVOBASE[[23]]$Means[c(1:3,6)]
m2 = EVOBASE[[24]]$Means[c(1:3,6)]
names(m1)==names(m2)

g1=meanStdG(g1, m1)*100
c(names(m1)==colnames(g1))

delta = m1-m2
delta = delta/sqrt(sum(delta^2)) #Unit length

evolvabilityBeta(g1, delta)$e
evolvabilityBeta(g1, delta)$c
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[3]*100

points(7, evolvabilityBeta(g1, delta)$e, pch=16)
points(7, evolvabilityMeans(g1)[1])
points(7, evolvabilityMeans(g1)[2])
points(7, evolvabilityMeans(g1)[3])

180-acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector

g2 = droptraits(EVOBASE[[24]]$G)

m1 = EVOBASE[[23]]$Means[c(1:3,5:6)]
m2 = EVOBASE[[24]]$Means[c(1:3,5:6)]
names(m1)==names(m2)

g2=meanStdG(g2, m2)*100
c(names(m2)==colnames(g2))

delta = m1-m2
delta = delta/sqrt(sum(delta^2)) #Unit length

evolvabilityBeta(g2, delta)$e
evolvabilityBeta(g2, delta)$c
evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[3]*100

points(8, evolvabilityBeta(g2, delta)$e, pch=16)
points(8, evolvabilityMeans(g2)[1])
points(8, evolvabilityMeans(g2)[2])
points(8, evolvabilityMeans(g2)[3])

acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G2 vs. divergence vector

#### Ipomopsis ####
#gsp
g1 = droptraits(EVOBASE[[18]]$G)[1:3, 1:3]
g1 = meanStdG(g1, EVOBASE[[18]]$Means[1:3])*100

#dsp
m1 = POPBASE[[35]]$popmeans[1, 2:4]
m2 = POPBASE[[35]]$popmeans[2, 2:4]
c(names(m1)==colnames(g1))

delta = c(as.matrix(m1-m2))
delta = delta/sqrt(sum(delta^2)) #Unit length

evolvabilityBeta(g1, delta)$e
evolvabilityBeta(g1, delta)$c
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[3]*100

180-acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector

#### Lobelia ####
#gsp
g1 = droptraits(EVOBASE[[1]]$G)
g1 = meanStdG(g1, EVOBASE[[1]]$Means)*100
g2 = droptraits(EVOBASE[[2]]$G)
g2 = meanStdG(g2, EVOBASE[[2]]$Means)*100

#dsp
m1 = POPBASE[[7]]$popmeans[1,c(6,7,4,5,2,3)]
m2 = POPBASE[[7]]$popmeans[2,c(6,7,4,5,2,3)]
c(names(m1)==colnames(g1))

delta = c(as.matrix(m1-m2))
delta = delta/sqrt(sum(delta^2)) #Unit length

evolvabilityBeta(g1, delta)$e
evolvabilityBeta(g1, delta)$c
evolvabilityMeans(g1)

evolvabilityBeta(g2, delta)$e
evolvabilityBeta(g2, delta)$c
evolvabilityMeans(g2)

acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector
acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G2 vs. divergence vector


# All the studies ####
nG = NULL
nD = NULL
for(s in c(1:length(both_sp))){
  species = both_sp[s]  
  nG[s] = length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD[s] = length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
}
cbind(both_sp, nG, nD)

reslist = list()
for(s in 1:length(both_sp)){
  for(g in 1:nG[s]){
    for(d in 1:nD[s]){
      res = computeGD(species = both_sp[s], gmatrix = g, dmatrix = d)[c(1:2, 5:18)]
      reslist[length(reslist)+1] = as.data.frame(unlist(res)[c(1:16)])
    }
  }
}
length(reslist)
reslist[[20]]

resmat = as.data.frame(matrix(NA, ncol=16, nrow=length(reslist)))
for(i in 1:length(reslist)){
  resmat[i,1:2]=as.character(reslist[[i]][1:2])
  resmat[i,3:13]=round(as.numeric(as.character(reslist[[i]][3:13])),2)
  resmat[i,14:16]=as.numeric(as.character(reslist[[i]][14:16]))
}

colnames(resmat)=c("D", "G", "npop", "theta", "betaG", "r2G", "betaD", "r2D","betaT","r2T","i_mean","e_mean","d_mean","nBetaG","nBetaD","nBetaT")

for(i in 1:nrow(resmat)){
  if(resmat$theta[i]>=90){
    resmat$theta[i] = 180 - resmat$theta[i]
  }
}

names(resmat)
#resmat = resmat[,c(2,1,16,3, 12,13,4,5,6)]
head(resmat)
View(resmat)

plot(resmat$theta, resmat$betaG, pch=16)
plot(resmat$theta, resmat$r2G, pch=16)

resmat = resmat[resmat$npop>2,]
resmat = resmat[resmat$nBetaG>2,]
resmat = resmat[resmat$nBetaD>2,]
resmat = resmat[resmat$nBetaT>2,]

mean(resmat$betaG)
hist(resmat$r2G)
hist(resmat$betaG)

plot(resmat$betaT, resmat$betaG)
plot(resmat$r2T, resmat$r2G)

plot(resmat$theta, resmat$betaG)

plot(resmat$npop, resmat$betaG)
plot(resmat$npop, resmat$r2G, cex=resmat$nBetaG*.5)

plot(resmat$i_mean, resmat$r2G, cex=resmat$nBetaG*.5, xlim=c(0,1))

#### Error-corrected D original code ####

drop = which(is.na(rowSums(means[-1])))
if(length(drop)>0){
  means=means[-drop]
  eV=ev[-drop]
}

#Set prior
n = ncol(means)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1))

means[,1:ncol(means)] = apply(means[,1:ncol(means)], 2, function(x) x*10)
mev = melt(eV)$value*100
data = means
vars = paste0(colnames(means), collapse=", ")
vars
vars = paste0("c(",vars,") ~-1+trait")

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

mod<-MCMCglmm(as.formula(noquote(vars)),
              rcov = ~us(trait):units,
              mev = mev,
              data = means, 
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

modD = matrix(apply(mod$VCV, 2, median)[2:(1+(ncol(means))^2)], nrow=n)
colnames(modD) = rownames(modD) = colnames(means)
modD = meanStdG(modD, colMeans(means))
round(modD, 3)
