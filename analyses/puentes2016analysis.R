##############################
#### - Example analyses - ####
##############################
rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)


#### Data from Puentes et al. 2016
#### Raw phenotypic data + breeding design for multiple populations
Gdat = read.csv("Z:/data/Evoldiv/data/puentes/puentes_etal_phen_multivariate_G_matrices.csv")

head(Gdat)
popmeans=apply(Gdat[,5:12], 2, function(x) tapply(x, Gdat$pop, mean, na.rm = T))
popmeans=popmeans[,c(2,3,5,6)]
popmeans

vars=apply(Gdat[,5:12], 2, function(x) tapply(x, Gdat$pop, var, na.rm=T))
vars=vars[,c(2,3,5,6)]
vars

Dmat=cov(log(popmeans))
Dmat



# G matrix for STUC population 
stuc = subset(Gdat, pop=="STUC")
stuc = droplevels(stuc)

stuc$cID=1:nrow(stuc)
stuc$animal=paste(as.character(stuc$sire), as.character(stuc$dam), stuc$cID, sep="_")

stuc$rosette.size.cm=sqrt(stuc$rosette.size.cm2)

#Mean-scale and multiply by 10
stuc[,c(6,7,9,15)]=apply(stuc[,c(6,7,9,15)], 2, function(x) 10*x/mean(x, na.rm=T))

ped=subset(stuc, select = c("animal", "dam", "sire"))
parentped <- cbind(animal=unique(c(as.character(stuc$sire), as.character(stuc$dam))), dam=NA, sire=NA)
ped<-rbind(parentped, ped)
#ped <- MCMCglmm:: prunePed(ped, stuc$animal) 
head(ped)

invA = inverseA(ped)$Ainv

n = 4
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 1000
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(petal.width.mm, petal.length.mm, flowers, rosette.size.cm) ~ -1+trait,
                random = ~us(trait):animal,
                rcov = ~us(trait):units,
                data = stuc, ginverse = list(animal = invA),
                family = rep("gaussian", n), prior = prior, 
                nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="Z:/data/Evoldiv/analyses/puentes2016/Gmat_STUC.RData")
#150k in 23.3 min

plot(mod$VCV[,12])
summary(mod$VCV)

gmat=matrix(apply(mod$VCV, 2, mean)[1:16],nrow=4)
colnames(gmat)=rownames(gmat)=c("petal.width.mm", "petal.length.mm", "flowers", "rosette.size.cm")
gmat

