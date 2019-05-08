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

Gmat=matrix(NA,nrow=6,ncol=6)
for(i in 1:6){
  for(j in 1:6){
    Gmat[i,j]=Gmat[j,i]=tmp[i,j]+tmp[j,i]
}}
diag(Gmat)=diag(Gmat)/2

Gmat

means=apply(Mik[,-1], 2, mean)
means
