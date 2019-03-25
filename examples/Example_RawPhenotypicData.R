##############################
#### - Example analyses - ####
##############################
rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)

#### Phenotypic data from Bartkowska et al. 2018
#### Raw phenotypic data, allows computing P for each population + D

indat=read.csv("data/eckert/Hierstruct_floral_morphology.csv")
pops=unique(indat$OutcropSystem)

#Compute P matrices and mean-scaled P-matrices per population
Plist=list()
msPlist=list()

for(i in 1:length(pops)){
  pop=pops[i]
  df=subset(indat[indat$OutcropSystem==pop,], 
            select = c("Herkogamy.mm","PistilLength.mm","SepalLength.mm","SpurLength.mm"))
  df=na.omit(df)
  P=cov(df)
  means=apply(df,2,mean,na.rm=T)
  msP=meanStdG(P,means)
  Plist[[i]]=P
  msPlist[[i]]=msP
  }

meanP=apply(simplify2array(Plist),1:2,mean)
meanP

mean_msP=apply(simplify2array(msPlist),1:2,mean)
mean_msP

popmeans=ddply(indat,.(OutcropSystem),summarize,
               z1=mean(Herkogamy.mm,na.rm=T),
               z2=mean(PistilLength.mm,na.rm=T),
               z3=mean(SepalLength.mm,na.rm=T),
               z4=mean(SpurLength.mm,na.rm=T))

popse2=ddply(indat,.(OutcropSystem),summarize,
            z1se=var(Herkogamy.mm,na.rm=T)/sum(Herkogamy.mm>-1,na.rm=T),
            z2se=var(PistilLength.mm,na.rm=T)/sum(PistilLength.mm>-1,na.rm=T),
            z3se=var(SepalLength.mm,na.rm=T)/sum(SepalLength.mm>-1,na.rm=T),
            z4se=var(SpurLength.mm,na.rm=T)/sum(SpurLength.mm>-1,na.rm=T))
popse2


D=cov(popmeans[,2:5])

#Compute error variance-covariance from P and n
meanN=mean(tapply(indat$OutcropSystem, indat$OutcropSystem, length))
De=meanP/meanN

#Compute and scale error-corrected D matrix
Dc=D-De
Dmatscaled=meanStdG(Dc,apply(popmeans[,2:5],2,mean))

#Gmax
ev1=eigen(mean_msP)$vectors[,1]
gmax=evolvabilityBeta(mean_msP, Beta = ev1)$e
dmax=evolvabilityBeta(Dmatscaled, Beta = ev1)$e

#Gmin
ev4=eigen(mean_msP)$vectors[,4]
gmin=evolvabilityBeta(mean_msP, Beta = ev4)$e
dmin=evolvabilityBeta(Dmatscaled, Beta = ev4)$e

#Random betas
betas=randomBeta(1000,4)
ebeta=evolvabilityBeta(mean_msP, betas)$e
dbeta=evolvabilityBeta(Dmatscaled, betas)$e

#Plot
plot(ebeta,dbeta,col="grey")
points(diag(mean_msP),diag(Dmatscaled),pch=16)
points(gmax,dmax,col="red", pch=16)
points(gmin,dmin,col="blue", pch=16)




#Mixed-model####
head(indat)

n = 4
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 200
burnin = samples*thin*.5
nitt = (samples*thin)+burnin


popmeans[,2:5]=apply(popmeans[,2:5],2,function(x)x*100)
mev=c(popse2$z1se, popse2$z2se, popse2$z3se, popse2$z4se)*10000

mod<-MCMCglmm(c(z1, z2, z3, z4) ~ -1+trait,
              random = ~us(trait):OutcropSystem,
              rcov = ~us(trait):units,
              mev = mev,
              data = popmeans, 
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

summary(mod)

modD=matrix(apply(mod$VCV, 2, mean)[1:16]/10000,nrow=4)
colnames(modD) = rownames(modD) = c("Herkogamy.mm","PistilLength.mm","SepalLength.mm","SpurLength.mm")
modD

cov(popmeans[,-1]/100)
Dc
