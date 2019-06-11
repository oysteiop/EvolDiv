##############################
#### - Example analyses - ####
##############################
rm(list=ls())
library(evolvability)
library(plyr)
library(reshape2)
library(MCMCglmm)
library(lme4)
library(kinship2)



#### Nigella data from Stefan Andersson
#### Mid-parent and offspring values from 2 populations
rm(list=ls())
list.files(path="./data/andersson/Nigella")
dat = read.csv2("./data/andersson/Nigella/Nigella_edited.csv", dec=".")
head(dat)
str(dat)

#####################################
#### - Estimating the G matrix - ####
#####################################

Mik=dat[đat$Population=="Mikonos",]
summary(lm(Mik[,8]~Mik[,2], na=na.exclude, data=Mik))$coef[2,1]

tmp=matrix(NA,nrow=6,ncol=6)
for(i in 1:6){
  for(j in 1:6){
    tmp[i,j]=cov(Mik[,i+1],Mik[,j+7])
}}
tmp

gmat=matrix(NA,nrow=6,ncol=6)
for(i in 1:6){
  for(j in 1:6){
    gmat[i,j]=gmat[j,i]=tmp[i,j]+tmp[j,i]
}}
diag(gmat)=diag(gmat)/2

colnames(gmat)=rownames(gmat)=c("Plant height", "Leaf distance", "Leaf length", "Sepal length", "Nectary length", "Anther length")
gmat

means=apply(Mik[,-1], 2, mean)[1:6]
means

meanStdG(gmat, means)*100

#####################################
#### - Estimating the D matrix - ####
#####################################

rm(list=ls())
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

#Compute eigenvectors etc.
first_ev=eigen(gmat)$vectors[,1]
gmax=evolvabilityBeta(gmat, Beta = first_ev)$e
dmax=evolvabilityBeta(dmat, Beta = first_ev)$e

last_ev=eigen(gmat)$vectors[,nrow(gmat)]
gmin=evolvabilityBeta(gmat, Beta = last_ev)$e
dmin=evolvabilityBeta(dmat, Beta = last_ev)$e

betas=randomBeta(1000,nrow(gmat))
ebeta=evolvabilityBeta(gmat, betas)$e
dbeta=evolvabilityBeta(dmat, betas)$e

plot(log10(ebeta),log10(dbeta),col="grey", las = 1,
     xlab="Evolvability",
     ylab="Population divergence",
     xlim=c(-1,2))
points(log10(diag(gmat)),log10(diag(dmat)),pch=16)
points(log10(gmax),log10(dmax),col="red", pch=16)
points(log10(gmin),log10(dmin),col="blue", pch=16)

