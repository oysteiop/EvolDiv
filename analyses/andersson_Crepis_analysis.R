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



#### Crepis data from Stefan Andersson
#### Among-pop and within-pop full sibs + selfed sibs
list.files(path="./data/andersson")
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
head(dat)
str(dat)

#####################################
#### - Estimating the G matrix - ####
#####################################
dat=dat[dat$TYPE=="OUTX",]
dat$animal=1:nrow(dat)
head(dat)

#Build the pedigree
pedigree=data.frame(as.character(dat$animal))
pedigree$dam=paste0(dat$IDNO,"d")
pedigree$sire=paste0(dat$IDNO,"s")
names(pedigree)=c("animal","dam","sire")

parentped=data.frame(animal=c(unique(pedigree$dam), unique(pedigree$sire)))
parentped$dam=rep(NA,nrow(parentped))
parentped$sire=rep(NA,nrow(parentped))
names(parentped)=c("animal","dam","sire")

pedigree=rbind(parentped, pedigree)
head(pedigree)

#Mean-scale and multiply by 10
#names(dat)[3]="animal"
dat[,c(4:8)]=apply(dat[,c(4:8)], 2, function(x) 10*x/mean(x, na.rm=T))
head(dat)

#Simple analysis
m=lmer(LEN~1+(1|IDNO),data=dat)
summary(m)
.568*2
(20.29*2)/63.13^2*100

m=lmer(MAX~1+(1|IDNO),data=dat)
summary(m)
.6205*2

m=lmer(MIN~1+(1|IDNO),data=dat)
summary(m)
.6205*2


#Run the MCMCglmm analysis
#Aped <- 2 * kinship2::kinship(pedigree[, 1], pedigree[,2], pedigree[, 3])

invA = inverseA(pedigree)$Ainv

#Five traits
n = 5
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(LEN,TIP,MAX,MIN,TEETH) ~ -1+trait,
              random = ~us(trait):animal,
              rcov = ~us(trait):units,
              data = dat, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/andersson_crepis/Gmat75k.RData")

#Check convergence
summary(mod$VCV)
plot(mod$VCV[,19])

#####################################
#### - Estimating the D matrix - ####
#####################################
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")
dat=dat[dat$TYPE=="POP",]
head(dat)

dat[,c(4:8)]=apply(dat[,c(4:8)], 2, function(x) log(x))
head(dat)

#Simple analysis
m=lmer(LEN~1+(1|IDENTITY),data=dat)
summary(m)

m=lmer(MAX~1+(1|IDNO),data=dat)
summary(m)

m=lmer(MIN~1+(1|IDNO),data=dat)
summary(m)


#Five traits
n = 5
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*5
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(LEN,TIP,MAX,MIN,TEETH) ~ -1+trait,
              random = ~us(trait):IDENTITY,
              rcov = ~us(trait):units,
              data = dat,
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/andersson_crepis/Dmat75k.RData")

#Check convergence
summary(mod$VCV)
plot(mod$VCV[,19])

################################
#### -Divergence analysis - ####
################################
dat = read.csv2("./data/andersson/Crepis_leaf_data.csv", dec=".")

#The G matrix
load(file="./analyses/andersson_crepis/Gmat75k.RData")
n=5
gmat=matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat)=rownames(gmat)=colnames(dat)[4:8]
gmat

#The D matrix
load(file="./analyses/andersson_crepis/Dmat75k.RData")

dmat=matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat)=rownames(dmat)=colnames(dat)[4:8]
dmat=dmat*100
dmat

#Compute eigenvectors etc.
first_ev=eigen(gmat)$vectors[,1]
gmax=evolvabilityBeta(gmat, Beta = first_ev)$e
dmax=evolvabilityBeta(gmat, Beta = first_ev)$e

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

