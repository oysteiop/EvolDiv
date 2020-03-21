###########################################################
#### Analyses of G and D matrices #####
###########################################################

library(reshape2)
library(MCMCglmm)
library(evolvability)
source("code/computeGD.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")

#Species present in both databases
gsp = unlist(lapply(EVOBASE, function(x) x$Species))
dsp = unlist(lapply(POPBASE, function(x) x$Species))
both_sp = unique(gsp[which(gsp %in% dsp)])
both_sp

POPBASE = POPBASE[-c(14:15)] #Dropping second M. guttatus study

out = computeGD(species = both_sp[13], gmatrix = 2, dmatrix = 1)
out

#### Estimate error-corrected D matrix ####
colnames(out$D)
means = POPBASE[[30]]$popmeans
means = means[,c(1,5,4,6,3,2,7,8)]
eV = POPBASE[[30]]$eV
eV = eV[,c(1,5,4,6,3,2,7,8)]

colnames(out$D)==names(means)[-1]

drop = which(is.na(rowSums(means[-1])))
if(length(drop)>0){
  means=means[-drop]
  eV=ev[-drop]
}

#Set prior
n = ncol(means)-1
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1))

means[,2:ncol(means)] = apply(means[,2:ncol(means)], 2, function(x) x*10)
mev = melt(eV[,-1])$value*100
data = means[,-1]
vars = paste0(colnames(means)[-1], collapse=", ")
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
#summary(mod)

modD = matrix(apply(mod$VCV, 2, median)[2:(1+(ncol(means)-1)^2)], nrow=n)
colnames(modD) = rownames(modD) = colnames(means)[-1]
modD = meanStdG(modD, colMeans(means[,-1]))
round(modD, 3)

eigen(modD)
plot(modD, out$D)
lines(-1:1, -1:1)
cor(c(modD), c(out$D))

#out$D = out$D*100
out$D = modD*100
out$G = out$G*100

evolvabilityMeans(out$G)
evolvabilityMeans(out$D)

#Compute eigenvectors etc.
g_ev = eigen(out$G)$vectors
var_g_g = evolvabilityBeta(out$G, Beta = g_ev)$e
var_d_g = evolvabilityBeta(out$D, Beta = g_ev)$e

d_ev = eigen(out$D)$vectors
var_g_d = evolvabilityBeta(out$G, Beta = d_ev)$e
var_d_d = evolvabilityBeta(out$D, Beta = d_ev)$e

#Compute summary stats
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

#Plot
x11(width=5, height=5)
xmin=-1.5
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin=-6
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Brassica cretica: Go", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(out$G)), log10(diag(out$D)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))


# All the studies
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


#####################################################
#### Evolvability along the vector of divergence ####
#####################################################

# Mimulus guttatus
g1 = droptraits(EVOBASE[[3]]$G)[1:2, 1:2]
g2 = droptraits(EVOBASE[[4]]$G)[1:2, 1:2]
g1 = droptraits(EVOBASE[[3]]$G)
g2 = droptraits(EVOBASE[[4]]$G)


m1 = EVOBASE[[3]]$Means[c(2,4:9)][1:2]
m2 = EVOBASE[[4]]$Means[c(2,4:9)][1:2]
m1 = EVOBASE[[3]]$Means[c(2,4:9)]
m2 = EVOBASE[[4]]$Means[c(2,4:9)]

g1=meanStdG(g1, m1)
g2=meanStdG(g2, m2)

delta = m1-m2
delta = delta/sqrt(sum(delta^2))

evolvabilityBeta(g1, delta)$e
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[3]*100

evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[3]*100

acos(t(eigen(g1)$vectors[,1]) %*% eigen(g2)$vectors[,1])*(180/pi) #Angle between the 2 G matrices
180-acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector
180-acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G2 vs. divergence vector

# Ellipse plot
ctr    <- c(0,0)  # data centroid
angles <- seq(0, 2*pi, length.out=500) 
A = g2
eigVal  <- eigen(A)$values[1:2]
eigVec  <- eigen(A)$vectors[1:2, 1:2]
eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
ellRot  <- eigVec %*% t(ellBase)                                          # rotated ellipse

plot((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, 
     xlab="Trait 1", ylab="Trait 2", xlim=c(-.06, .06), ylim=c(-.06,.06),
     type="l", lwd=1, col="lightblue", las=1)
polygon((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, col="lightblue")
matlines(xMat[,1:2], yMat[,1:2], lty=1, lwd=1, col="black")

ctr    <- div[1:2]  # data centroid
angles <- seq(0, 2*pi, length.out=500) 
A = g2
eigVal  <- eigen(A)$values[1:2]
eigVec  <- eigen(A)$vectors[1:2, 1:2]
eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
ellRot  <- eigVec %*% t(ellBase)                                          # rotated ellipse

polygon((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, col="lightblue")
matlines(xMat[,1:2], yMat[,1:2], lty=1, lwd=1, col="black")

delta
div = colSums(delta * t(eigen(g1)$vectors))
arrows(0,0, div[1], div[2], length=.1)



library(ellipse)
lines(ellipse(g1, which=c(1,7), level=.35), lwd=3)

# Nigella degenii
g1 = droptraits(EVOBASE[[3]]$G)[1:2, 1:2]
g2 = droptraits(EVOBASE[[4]]$G)[1:2, 1:2]
g1 = droptraits(EVOBASE[[23]]$G)
g2 = droptraits(EVOBASE[[24]]$G)


m1 = EVOBASE[[23]]$Means[c(1:3,5:6)]
m2 = EVOBASE[[24]]$Means[c(1:3,5:6)]

g1=meanStdG(g1, m1)
g2=meanStdG(g2, m2)


delta = m1-m2
delta = delta/sqrt(sum(delta^2))

evolvabilityBeta(g1, delta)$e
evolvabilityMeans(g1)
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[3]*100

evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[3]*100

180-acos(t(eigen(g1)$vectors[,1]) %*% eigen(g2)$vectors[,1])*(180/pi) #Angle between the 2 G matrices
acos(t(eigen(g1)$vectors[,1]) %*% delta)*(180/pi) #G1 vs. divergence vector
180-acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi) #G2 vs. divergence vector
