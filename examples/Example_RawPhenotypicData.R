##############################
#### - Example analyses - ####
##############################

library(evolvability)
library(plyr)


#### Phenotypic data from Bartkowska et al. 2018
#### Raw phenotypic data, allows computing P for each population + D

indat=read.csv("data/eckert/Hierstruct_floral_morphology.csv")
names(indat)

pops=unique(indat$OutcropSystem)

outlist=list()
for(i in 1:length(pops)){
  pop=pops[i]
  df=subset(indat[indat$OutcropSystem==pop,], 
            select = c("Herkogamy.mm","PistilLength.mm","SepalLength.mm","SpurLength.mm"))
  df=na.omit(df)
  P=cov(df)
  means=apply(df,2,mean,na.rm=T)
  msP=meanStdG(P,means)
  outlist[[i]]=msP
}
outlist

meanP=apply(simplify2array(outlist),1:2,mean)
meanP

popmeans=ddply(indat,.(OutcropSystem),summarize,
               z1=mean(Herkogamy.mm,na.rm=T),
               z2=mean(PistilLength.mm,na.rm=T),
               z3=mean(SepalLength.mm,na.rm=T),
               z4=mean(SpurLength.mm,na.rm=T))

popse=ddply(indat,.(OutcropSystem),summarize,
            z1se=sd(Herkogamy.mm,na.rm=T)/sqrt(sum(Herkogamy.mm>-1,na.rm=T)),
            z2se=sd(PistilLength.mm,na.rm=T)/sqrt(sum(PistilLength.mm>-1,na.rm=T)),
            z3se=sd(SepalLength.mm,na.rm=T)/sqrt(sum(SepalLength.mm>-1,na.rm=T)),
            z4se=sd(SpurLength.mm,na.rm=T)/sqrt(sum(SpurLength.mm>-1,na.rm=T)))
popse

D=cov(popmeans[,2:5])
Dmatscaled=meanStdG(D,apply(popmeans[,2:5],2,mean))

#Gmax
ev1=eigen(meanP)$vectors[,1]
gmax=evolvabilityBeta(meanP, Beta = ev1)$e
dmax=evolvabilityBeta(Dmatscaled, Beta = ev1)$e

#Gmin
ev4=eigen(meanP)$vectors[,4]
gmin=evolvabilityBeta(meanP, Beta = ev4)$e
dmin=evolvabilityBeta(Dmatscaled, Beta = ev4)$e

#Random betas
betas=randomBeta(1000,4)
ebeta=evolvabilityBeta(meanP, betas)$e
dbeta=evolvabilityBeta(Dmatscaled, betas)$e

#Plot
plot(ebeta,dbeta,col="grey")
points(diag(meanP),diag(Dmatscaled),pch=16)
points(gmax,dmax,col="red", pch=16)
points(gmin,dmin,col="blue", pch=16)


