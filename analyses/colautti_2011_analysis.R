###############################################
#### - Data from Colautti & Barrett 2011 - ####
###############################################

# Lythrum salicaria

rm(list=ls())
library(evolvability)
library(plyr)
library(reshape2)
library(MCMCglmm)

#Gdat = read.table("data/colautti/AllDataFixed.txt", dec=".", sep=",", header=T)
Gdat = read.table("data/colautti/AllDataFixed2.txt", header=T)
names(Gdat)

popmeans = as.data.frame(apply(Gdat[,c(11:23, 25:38)], 2, function(x) tapply(x, Gdat$Pop, mean, na.rm = T)))
head(popmeans)

popmeans = subset(popmeans, select=c("TLeafArea", "THeight", "Height2wk", "Height4wk",
                                     "FDays", "FStemWidth", "FVeg",
                                     "FInf", "HVeg","HInf","HVegW","HInfW"))
head(t(popmeans))
write.csv2(melt(t(popmeans)), file = "data/colautti/popmeans.csv", row.names=T)

popvars = as.data.frame(apply(Gdat[,c(11:23, 25:38)], 2, function(x) tapply(x, Gdat$Pop, var, na.rm = T)))
head(popvars)

popvars = subset(popvars, select=c("TLeafArea", "THeight", "Height2wk", "Height4wk",
                                   "FDays", "FStemWidth", "FVeg",
                                   "FInf", "HVeg","HInf","HVegW","HInfW"))
head(t(popvars))
write.csv2(melt(t(popvars)), file = "data/colautti/popvars.csv", row.names=T)

popn = as.data.frame(apply(Gdat[,c(11:23, 25:38)], 2, function(x) tapply(x>0, Gdat$Pop, sum, na.rm = T)))
head(popn)

popn = subset(popn, select=c("TLeafArea", "THeight", "Height2wk", "Height4wk",
                             "FDays", "FStemWidth", "FVeg",
                             "FInf", "HVeg","HInf","HVegW","HInfW"))
head(t(popn))
write.csv2(melt(t(popn)), file = "data/colautti/popn.csv", row.names=T)

dmat = cov(log(popmeans))*100
dmat

means = apply(Gdat[,c(11:23, 25:38)], 2, mean, na.rm=T)
means
cbind(means[c(23,1,5,7,9,14,11,12,16,17,19,20)])

vars = apply(Gdat[,c(11:23, 25:38)], 2, var, na.rm=T)
vars
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
plot(ev, c(diag(dmat), var(popmeans$logVeg, na.rm=T)*100, var(popmeans$LogInf, na.rm=T)*100))


#### G matrix for one population #### 
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

# All
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

rownames(emat)=pops
colnames(emat)=traits

evals=colMeans(emat)

evals

#####################################
######## - Analyse G vs. D - ########
#####################################

#Univariate
plot(evals, diag(dmat))
plot(log10(evals), log10(diag(dmat)))

