##############################
#### - Example analyses - ####
##############################
rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)


#### Data from Yu et al. 2011
#### Breeding design with within and among-pop crosses, 3 populations

list.files(path="./data/delph")
Gdat = read.table("./data/delph/traitdata.txt", header=T)
crossID=read.table("./data/delph/crossID.txt", header=T)

#Format the trait data
calyx_w3_F=rep(NA, nrow(Gdat))
calyx_w3_M=rep(NA, nrow(Gdat))
calyx_w4_F=rep(NA, nrow(Gdat))
calyx_w4_M=rep(NA, nrow(Gdat))
calyx_w5_F=rep(NA, nrow(Gdat))
calyx_w5_M=rep(NA, nrow(Gdat))

for(i in 1:nrow(Gdat)){
  if(Gdat$sex[i]=="F"){
    calyx_w3_F[i]=Gdat$calyx_w3[i]
    calyx_w4_F[i]=Gdat$calyx_w4[i]
    calyx_w5_F[i]=Gdat$calyx_w5[i]
  }
  else{
    calyx_w3_M[i]=Gdat$calyx_w3[i]
    calyx_w4_M[i]=Gdat$calyx_w4[i]
    calyx_w5_M[i]=Gdat$calyx_w5[i]
  }}

trdata=data.frame(calyx_w3_F,calyx_w3_M,calyx_w4_F,calyx_w4_M,calyx_w5_F,calyx_w5_M)
trdata$mean_calyx_F=rowMeans(subset(trdata, select=c("calyx_w3_F", "calyx_w4_F","calyx_w5_F")), na.rm=T)
trdata$mean_calyx_M=rowMeans(subset(trdata, select=c("calyx_w3_M", "calyx_w4_M","calyx_w5_M")), na.rm=T)

data=data.frame(Gdat[,1:3],trdata)
names(data)[1]="Cross"
data=data[-1020,] #Remove duplicated ID

crossID=crossID[,1:3]
data=merge(data, crossID, by="Cross", all.x=T)
data$animal=paste(data$Cross, data$individual, sep="-")

head(data)

#Format the pedigree data
pedigree=subset(data, select=c("animal","Female","Male"))
names(pedigree)=c("animal","dam","sire")
parentped=data.frame(animal=unique(c(as.character(pedigree$dam), as.character(pedigree$sire))), dam=NA,sire=NA)
pedigree=rbind(parentped, pedigree)
pedigree$animal=factor(pedigree$animal)
pedigree$dam=factor(pedigree$dam)
pedigree$sire=factor(pedigree$sire)
head(pedigree,30)

#Mean-scale and multiply by 10
head(data)
data[,c(4:11)]=apply(data[,c(4:11)], 2, function(x) 10*x/mean(x, na.rm=T))

#Run the MCMCglmm analysis
invA = inverseA(pedigree)$Ainv

#Two traits
n = 2
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 10
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(mean_calyx_F, mean_calyx_M) ~ -1+trait,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = data, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/delph_Silene/Gmat15k_2traits.RData")

#Six traits
n = 2
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 10
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(calyx_w3_F, calyx_w4_F, calyx_w5_F, calyx_w3_M, calyx_w4_M, calyx_w5_M) ~ -1+trait,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = data, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/delph_Silene/Gmat15k.RData")


load(file="./analyses/delph_Silene/Gmat15k.RData")

#Check convergence
summary(mod$VCV)
plot(mod$VCV[,1])

#Compile G matrix
gmat=matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat)=rownames(gmat)=c("calyx_w3_F", "calyx_w4_F", "calyx_w5_F", "calyx_w3_M", "calyx_w4_M", "calyx_w5_M")
colnames(gmat)=rownames(gmat)=c("calyx_F", "calyx_M")
gmat

evolvabilityBeta(gmat[1:3,1:3], Beta=c(1,0,0))
evolvabilityBeta(gmat[4:6,4:6], Beta=c(1,0,0))

# - D matrix - ####
head(data)
unique(data$Female)
data$dampop=substr(data$Female,1,1)
data$sirepop=substr(data$Male,1,1)

reddat=data[which(data$dampop==data$sirepop),]
head(reddat)
popmeans=ddply(reddat,.(dampop, sirepop), summarize,
               calyx_w3_F=mean(calyx_w3_F, na.rm=T),
               calyx_w4_F=mean(calyx_w4_F, na.rm=T),
               calyx_w5_F=mean(calyx_w5_F, na.rm=T),
               calyx_w3_M=mean(calyx_w3_M, na.rm=T),
               calyx_w4_M=mean(calyx_w4_M, na.rm=T),
               calyx_w5_M=mean(calyx_w5_M, na.rm=T))

popmeans

Dmat=cov(log(popmeans[,3:8]))*100
Dmat

plot(diag(gmat), diag(Dmat))

gmat=gmat[c(1,4),c(1,4)]
Dmat=Dmat[c(1,4),c(1,4)]

eigen(Dmat)
#Gmax
ev1=eigen(gmat)$vectors[,1]
gmax=evolvabilityBeta(gmat, Beta = ev1)$e
dmax=evolvabilityBeta(Dmat, Beta = ev1)$e

#Gmin
ev4=eigen(gmat)$vectors[,nrow(gmat)]
gmin=evolvabilityBeta(gmat, Beta = ev4)$e
dmin=evolvabilityBeta(Dmat, Beta = ev4)$e

#Random betas
betas=randomBeta(1000,nrow(gmat))
ebeta=evolvabilityBeta(gmat, betas)$e
dbeta=evolvabilityBeta(Dmat, betas)$e

#Plot
plot(ebeta,dbeta,col="grey")
points(diag(gmat),diag(Dmat),pch=16)
points(gmax,dmax,col="red", pch=16)
points(gmin,dmin,col="blue", pch=16)
