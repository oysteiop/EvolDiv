###############################################
#### - Univariate analysis of D database - ####
###############################################
rm(list=ls())

library(plyr)
library(reshape2)
library(lme4)
library(MCMCglmm)
#library(devtools)
#library(withr)
#with_libpaths(new="C:/Program Files/R/R-3.5.0/library",install_github("GHBolstad/evolvability"))
#.libPaths("C:/Program Files/R/R-3.5.0/library")
library(evolvability)

ddat = read.table("data/dmatdata.txt", header=T)
ddat$ID = paste(ddat$reference,ddat$species,ddat$environment, sep="_")
maxdists = read.csv(file="data/maxdists.csv")

#Scaling of trait variances with the mean
#ddat=ddat[ddat$trait!="organ_size",]
#plot(log10(ddat$mean), log10(ddat$sd))
#lines(-10:10, -10:10)

#List of studies
studies=sort(unique(ddat$ID))
studies
splist=tapply(as.character(ddat$species), ddat$ID, function(x) x[1])

#Compute among-pop variance for each trait in each study
outlist=list()
for(s in 1:length(studies)){
  red=ddat[ddat$ID==studies[[s]],]
  red$trait=factor(red$trait)
  df=data.frame(study_ID=studies[[s]], 
                species=splist[[s]], 
                trait=sort(unique(red$trait)),
                environment=sort(unique(red$environment)),
                npop=tapply(red$mean>-Inf, red$trait, sum, na.rm=T),
                d=cbind(tapply(log(red$mean), red$trait, var, na.rm=T)),
                mean=cbind(tapply(red$mean, red$trait, mean, na.rm=T)))
  df=na.omit(df)
  outlist[[s]]=df
}
ddf=rbind.fill(outlist)
head(ddf,5)

# Sampling variance of a variance (From Lynch & Walsh 1998 p. 815)
sv = (2*(ddf$d^2))/(ddf$npop+2)
ddf$d_se = sqrt(sv)

#Combine with data from Evolvability database
edat=read.table("data/evolvabilitydatabase2019.txt", header=T)

edat$species_measurement=paste0(edat$species,"_",edat$measurement)
ddf$species_trait=paste0(ddf$species,"_",ddf$trait)

tg1=NULL
for(i in 1:nrow(ddf)){
  tg1[i]=as.character(edat$traitgroup1)[which(as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])[1]]
}

tg2=NULL
for(i in 1:nrow(ddf)){
  tg2[i]=as.character(edat$traitgroup2)[which(as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])[1]]
}

dimension=NULL
for(i in 1:nrow(ddf)){
  dimension[i]=as.character(edat$dimension)[which(as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])[1]]
}

ms=NULL
for(i in 1:nrow(ddf)){
  ms[i]=as.character(edat$matingsystem)[which(as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])[1]]
}

maxdist=NULL
for(i in 1:nrow(ddf)){
  maxdist[i]=as.numeric(maxdists[which(maxdists[,2]==as.character(ddf$study_ID)[i])[1],3])
}

evals=NULL
for(i in 1:nrow(ddf)){
w=which(as.character(edat$species)==as.character(ddf$species)[i] & as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])
evals[i]=mean(edat$evolvability[w],na.rm=T)
}

ddf$maxdist=maxdist
ddf$tg1=tg1
ddf$tg2=tg2
ddf$ms=ms
ddf$ms=factor(ddf$ms, levels=c("S","M","O"))
ddf$dimension=dimension
ddf$dimension=factor(ddf$dimension, levels=c("linear", "area", "mass_volume", "count", "ratio"))
ddf$evals=evals
head(ddf)

ddf=ddf[ddf$d>0,]
#ddf=ddf[ddf$maxdist>0,]
ddf=ddf[ddf$evals>0,]
ddf=na.omit(ddf)
head(ddf)

#Check
tg1=tg2=dimension=ms=evals=NULL

for(i in 1:nrow(ddf)){
  w=which(as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])
  sel=edat[w,]
  
  tg1[i]=as.character(unique(sel$traitgroup1))
  tg2[i]=as.character(unique(sel$traitgroup2))
  dimension[i]=as.character(unique(sel$dimension))
  ms[i]=as.character(unique(sel$matingsystem))
  evals[i]=mean(sel$evolvability)
}


#### Test off-setting of area and cubic measures ####
for(i in 1:nrow(ddf)){
  if(ddf$dimension[i]=="area"){
    ddf$d[i]=ddf$d[i]/4
    ddf$evals[i]=ddf$evals[i]/4
  }
  if(ddf$dimension[i]=="mass_volume"){
    ddf$d[i]=ddf$d[i]/9
    ddf$evals[i]=ddf$evals[i]/9
  }
}

# Informal meta-analysis

#Remove repeated D. scandens studies?
#ddf=ddf[ddf$study_ID!="Hansen_et_al._2003_Dalechampia_scandens_A_greenhouse",]
#ddf=ddf[ddf$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_field",]
#ddf=ddf[ddf$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_greenhouse",]

m=lmer(log(d)~ log(evals) + log(npop) + log(maxdist) + ms + tg1 + dimension + (1|species/study_ID), data=ddf)
summary(m)

# Formal meta-analysis using Almer_SE

SE=sqrt((ddf$d_se^2)/(ddf$d^2))
plot(ddf$npop, SE)
m=Almer_SE(log(d) ~ log(evals) + log(npop) + log(maxdist) + ms + tg1 + dimension + (1|species/study_ID), 
           SE=SE, maxiter = 100, data=ddf)
summary(m)

#### Plotting evolvability vs. divergence ####

#All data
plot(log10(ddf$evals),log10(ddf$d*100),
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf$npop),
     col="black",
     xlim=c(-2.5, 2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

# Function to plot e vs. d with subset highlighted ####
plotSubset=function(category, subset, xlim=c(-2.5,2), ylim=c(-4,3), ...){
plot(log10(ddf$evals),log10(ddf$d*100),
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf$npop),
     col="grey", main=paste0(category, ": ",subset),
     xlim=xlim, ylim=ylim, xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3), 10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2), 10^c(-4,-3,-2,-1,0,1,2), las=1)

column=which(names(ddf)==category)
ddf2=ddf[ddf[,column]==subset,]
points(log10(ddf2$evals), log10(ddf2$d*100), col="black", cex=1*sqrt(ddf2$npop), ...)
}

x11()
par(mfrow=c(1,3))
plotSubset("ms", "S", lwd=2)
plotSubset("ms", "M", lwd=2)
plotSubset("ms", "O", lwd=2)

x11()
par(mfrow=c(2,3))
plotSubset("tg1", "floral")
plotSubset("tg1", "vegetative")
plotSubset("tg2", "flowersize")
plotSubset("tg2", "fit")
plotSubset("tg2", "display")
plotSubset("tg2", "reward")

x11()
par(mfrow=c(2,3))
plotSubset("dimension", "linear")
plotSubset("dimension", "area")
plotSubset("dimension", "mass_volume")
plotSubset("dimension", "count")
plotSubset("dimension", "ratio")

studies=unique(ddf$study_ID)

pdf("figs/univariate_G_D_plots2.pdf")
for(s in 1:length(studies)){
  study=studies[s]
  plotSubset("study_ID", study)
}
dev.off()

#### Formal meta-analysis using MCMCglmm ####

#Prepare data
names(ddf)
moddat=subset(ddf, select=c("d", "evals", "maxdist", "npop", "tg1", "dimension", "ms","study_ID", "species", "d_se", "mean"))
moddat=na.omit(moddat)
moddat$d=moddat$d*100
moddat$d_se=moddat$d_se*100

#moddat=moddat[moddat$dimension=="linear",]

#Set prior
prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002),
                                          G2=list(V=1, nu=0.002)))

prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002),
                                          G2=list(V=1, nu=0.002),
                                          G3=list(V = diag(1), fix = 1)))

#Compute mean-scaled measurement error variances
mev=moddat$d_se^2
moddat$SE=sqrt(mev)  

#To log-normal
moddat$test=(-2*log(moddat$d)) + log(mev+(moddat$d^2))
moddat$test=sqrt(moddat$test)

hist(moddat$test)
plot(moddat$npop, moddat$test)

#test2=mev/moddat$d^2

#Set sampling parameters
samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

#Sample MCMC
mod<-MCMCglmm(log(d) ~ log(evals) + log(maxdist) + log(npop) + ms + tg1 + dimension,
              random = ~study_ID + species + idh(test):units,
              rcov = ~units,
              #mev = mev,
              data = moddat, 
              family = "gaussian", prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

#Check convergence
x11()
plot(mod$VCV)

#Parameter estimates
summary(mod)


plot(log(moddat$maxdist), log(moddat$d))
plot(log(moddat$evals), log(moddat$d))
plot(log(moddat$npop), log(moddat$d))
