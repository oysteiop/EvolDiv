##############################
#### - Example analyses - ####
##############################
rm(list=ls())
library(evolvability)
library(plyr)


#### Data from Puentes et al. 2016
#### Raw phenotypic data + breeding design for multiple populations

Gdat = read.csv("data/puentes/puentes_etal_phen_multivariate_G_matrices.csv")
head(Gdat)
popmeans=apply(Gdat[,5:12], 2, function(x) tapply(x, Gdat$pop, mean, na.rm = T))
popmeans=popmeans[,c(2,3,5,6)]
popmeans

vars=apply(Gdat[,5:12], 2, function(x) tapply(x, Gdat$pop, var, na.rm=T))
vars=vars[,c(2,3,5,6)]
vars

Dmat=cov(log(popmeans))
Dmat


library(MCMCglmm)

# G matrix for STUC population 
genvar.stuc = subset(Gdat, pop=="STUC")
genvar.stuc = droplevels(genvar.stuc)

genvar.stuc$cID=1:nrow(genvar.stuc)
genvar.stuc$animal=paste(as.character(genvar.stuc$sire), as.character(genvar.stuc$dam), genvar.stuc$cID, sep="_")

ped=subset(genvar.stuc, select = c("animal", "dam", "sire"))
head(ped)

parentped <- cbind(animal=unique(c(as.character(genvar.stuc$sire), as.character(genvar.stuc$dam))), dam=NA, sire=NA)

ped<-rbind(parentped, ped)
#ped <- MCMCglmm:: prunePed(ped, genvar.stuc$animal) 

invA = inverseA(ped)$Ainv

nitt <- 30000
burnin <- 10000
thin <-20


n = 4
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

  mod<-MCMCglmm(c(petal.width.mm, petal.length.mm, flowers, rosette.size.cm2) ~ -1+trait,
                random= ~us(trait):animal,
                rcov=~us(trait):units,
                data=genvar.stuc, ginverse=list(animal=invA),
                family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)

gmat=matrix(apply(mod$VCV, 2, mean)[1:16],nrow=4)
colnames(gmat)=rownames(g2)=c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm2")
gmat
  
popmeans
scaledg=(diag(gmat)/popmeans[2,]^2)*100

plot(scaledg,log(diag(Dmat)))  
