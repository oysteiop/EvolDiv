################################################
#### - Nigella data from Stefan Andersson - ####
################################################

#### Mid-parent and offspring values from 2 populations

rm(list=ls())
library(evolvability)
library(plyr)
library(reshape2)
library(MCMCglmm)
library(lme4)
library(kinship2)

rm(list=ls())
list.files(path="./data/andersson/Nigella")
dat = read.csv2("./data/andersson/Nigella/Nigella_edited.csv", dec=".")
head(dat)
str(dat)

#####################################
#### - Estimating the G matrix - ####
#####################################

Mik=dat[dat$Population=="Mikonos",]

#Remove outlier
Mik=Mik[-which.max(Mik$mp_Leaf_length),]

#Log-transform
Mik[,2:13]=apply(Mik[,2:13], 2, log)

#x11()
#par(mfrow=c(3, 4))
#for(i in 2:13){hist(Mik[, i])}

#G-matrix from parent-offspring covariances
tmp = matrix(NA,nrow=6,ncol=6)
for(i in 1:6){
  for(j in 1:6){
    tmp[i,j] = cov(Mik[,i+1],Mik[,j+7])
}}
tmp

gmat = matrix(NA,nrow=6,ncol=6)
for(i in 1:6){
  for(j in 1:6){
    gmat[i,j] = gmat[j,i]=tmp[i,j]+tmp[j,i]
}}
#diag(gmat)=diag(gmat)/2
gmat=gmat/2

colnames(gmat)=rownames(gmat)=c("Plant height", "Leaf distance", "Leaf length", "Sepal length", "Nectary length", "Anther length")

gmat
cov2cor(gmat)

gmat=gmat*100
#means=apply(Mik[,-1], 2, mean)[1:6]
#means

#gmat=meanStdG(gmat, means)*100
#gmat
evolvabilityMeans(gmat)

#####################################
#### - Estimating the D matrix - ####
#####################################

list.files(path="./data/andersson/Nigella")
dat = read.csv2("./data/andersson/Nigella/Nigella_edited.csv", dec=".")
head(dat)
str(dat)

popmeans=ddply(dat, .(Population), summarize,
               Plant_height=mean(mp_Plant_height, na.rm=T),
               Leaf_distance=mean(mp_Leaf_distance, na.rm=T),
               Leaf_length=mean(mp_Leaf_length, na.rm=T),
               Sepal_length=mean(mp_Sepal_length, na.rm=T),
               Nectary_length=mean(mp_Nectary_length, na.rm=T),
               Anther_length=mean(mp_Anther_length, na.rm=T))
popmeans=popmeans[,-1]
dmat=cov(log(popmeans))

dmat=round(dmat, 10)
gmat=round(gmat, 10)

#Simple analysis
plot(log10(diag(gmat)), log10(diag(dmat)), pch=16)

#Compute eigenvectors etc.
g_ev=eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
#var_d_g = diag(t(g_ev) %*% dmat %*% g_ev)

d_ev=eigen(dmat)$vectors
var_g_d=evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d=evolvabilityBeta(dmat, Beta = d_ev)$e
#var_d_d = diag(t(d_ev) %*% dmat %*% d_ev)

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

#Plot
plot(var_g_g, var_d_g, xlim=c(-20,40), ylim=c(-.1,.15))
points(var_g_d, var_d_d, pch=16)
points(diag(gmat), diag(dmat), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), 
       pch=c(1,16, 16), col=c("black", "black", "blue"))


#On log10 scale
xmin=log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax=log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin=-3
ymax=log10(max(c(var_d_g, var_d_d), na.rm=T))

plot(log10(var_g_g), log10(var_d_g), xlim=c(-2,2), ylim=c(-12, 0))
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")

legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), 
       pch=c(1,16, 16), col=c("black", "black", "blue"))

#Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
