##############################
#### - Example analyses - ####
##############################
rm(list=ls())
library(evolvability)
library(plyr)


#### Data from Walter et al. 2018
#### Raw phenotypic data + breeding design for multiple populations


list.files()

#D-matrix
Ddat = read.csv("data/walter/Data_Exp1_Dmatrix.csv")
head(Ddat)
popmeans=apply(Ddat[,3:12], 2, function(x){tapply(x, Ddat$Population,mean,na.rm=T)})
popmeans
Dmat=cov(log(popmeans))
Dmat

#G-matrix
Gdat = read.csv("Data_Exp2_Gvariance.csv")
head(Gdat)
popmeans=apply(Gdat[,8:17], 2, function(x){tapply(x, Gdat$Type,mean,na.rm=T)})
popmeans

vars=apply(Gdat[,8:17], 2, function(x){tapply(x, Gdat$Type,var,na.rm=T)})
vars

(.294*4954.863)/(410.66^2)*100




