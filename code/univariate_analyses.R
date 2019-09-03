###############################################
#### - Univariate analysis of D database - ####
###############################################
rm(list=ls())
library(plyr)
library(reshape2)
library(evolvability)
library(lme4)

ddat=read.table("data/dmatdata.txt", header=T)
ddat$ID=paste(ddat$reference,ddat$species,ddat$environment, sep="_")

#List of studies
studies=sort(unique(ddat$ID))
studies
splist=tapply(as.character(ddat$species), ddat$ID, function(x) x[1])

#Compute among-pop variance for each trait in each study
outlist=list()
for(s in 1:length(studies)){
  red=ddat[ddat$ID==studies[[s]],]
  red$trait=factor(red$trait)
  #red$mean=log(red$mean)
  df=data.frame(study_ID=studies[[s]], 
                species=splist[[s]], 
                trait=sort(unique(red$trait)),
                npop=tapply(red$mean>-Inf, red$trait, sum, na.rm=T),
                d=cbind(tapply(log(red$mean), red$trait, var, na.rm=T)),
                mean=cbind(tapply(red$mean, red$trait, mean, na.rm=T)),
                dse2=cbind(tapply(red$se^2, red$trait, mean, na.rm=T)))
  df=na.omit(df)
  outlist[[s]]=df
}
ddf=rbind.fill(outlist)
head(ddf,5)

#Combine with data from Evolvability database
edat=read.table("data/evolvabilitydatabase2019.txt", header=T)

tg1=NULL
for(i in 1:nrow(ddf)){
  tg1[i]=as.character(edat$traitgroup1)[which(as.character(edat$measurement)==as.character(ddf$trait)[i])[1]]
}

dimension=NULL
for(i in 1:nrow(ddf)){
  dimension[i]=as.character(edat$dimension)[which(as.character(edat$measurement)==as.character(ddf$trait)[i])[1]]
}

evals=NULL
for(i in 1:nrow(ddf)){
w=which(as.character(edat$species)==as.character(ddf$species)[i] & as.character(edat$measurement)==as.character(ddf$trait)[i])
evals[i]=mean(edat$evolvability[w],na.rm=T)
}

ddf$tg1=tg1
ddf$dimension=dimension
ddf$evals=evals
head(ddf)

tapply(ddf$d, ddf$tg1, median, na.rm=T)*100
tapply(ddf$d>-Inf, ddf$tg1, sum, na.rm=T)
plot(as.factor(ddf$tg1), log10(ddf$d))

ddf$dimension=factor(ddf$dimension, levels=c("linear", "area", "mass_volume", "count", "ratio"))
levels(ddf$dimension)
tapply(ddf$d, ddf$dimension, median, na.rm=T)*100
tapply(ddf$d>-Inf, ddf$dimension, sum, na.rm=T)
plot(as.factor(ddf$dimension), log10(ddf$d))

#Informal meta-analysis
#Covariates: evolvability, n populations, max distance, traitgroup, dimension?, mating system?
#Mean within-pop se2 as measurement variance in formal meta-analysis

m=lmer(log(d)~npop + evals + tg1 + dimension+ (1|species/study_ID), data=ddf)
summary(m)

#Plot
ddf=ddf[ddf$d>0,]
ddf=ddf[ddf$evals>0,]
ddf=na.omit(ddf)

#ddf=ddf[ddf$species!="Dalechampia_scandens_A",]

plot(log10(ddf$evals),log10(ddf$d*100),
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf$npop),
     col=as.numeric(as.factor(ddf$tg1)),
     xlim=c(-2.5,2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

ddf2=ddf[ddf$dimension=="linear",]
plot(log10(ddf2$evals),log10(ddf2$d*100),
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf2$npop),
     col="blue",
     main="Linear traits",
     xlim=c(-2.5,2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

ddf2=ddf[ddf$dimension=="area",]
points(log10(ddf2$evals),log10(ddf2$d*100),col="red",cex=1*sqrt(ddf2$npop))
ddf2=ddf[ddf$dimension=="mass_volume",]
points(log10(ddf2$evals),log10(ddf2$d*100),col="darkgreen",cex=1*sqrt(ddf2$npop))
ddf2=ddf[ddf$dimension=="count",]
points(log10(ddf2$evals),log10(ddf2$d*100),col="red",cex=1*sqrt(ddf2$npop))
ddf2=ddf[ddf$dimension=="ratio",]
points(log10(ddf2$evals),log10(ddf2$d*100),col="black",cex=1*sqrt(ddf2$npop))

# Floral vs. vegetative traits
ddf2=ddf[ddf$tg1=="floral",]
plot(log10(ddf2$evals),log10(ddf2$d*100),
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf2$npop),
     col="blue",
     xlim=c(-2.5,2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

ddf2=ddf[ddf$tg1=="vegetative",]
points(log10(ddf2$evals),log10(ddf2$d*100),col="black",cex=1*sqrt(ddf2$npop))


sort(tapply(ddf$evals>-Inf, ddf$study_ID, sum, na.rm=T))

studies=unique(ddf$study_ID)

pdf("figs/univariate_G_D_plots.pdf")
for(s in 1:length(studies)){
  red=ddf[ddf$study_ID==studies[s],]
  plot(log(red$evals),log(red$d*100), pch=16, main=paste(studies[s]))
}
dev.off()


#Formal meta-analysis
#Can include e.g. trait groups, number of pops, max distance between pops,

#Prepare data
names(ddf)
moddat=subset(ddf, select=c("d", "evals", "tg1", "dimension", "study_ID", "species", "dse2", "mean"))
moddat=na.omit(moddat)
moddat$d=moddat$d*100

#Set prior
prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002),
                                          G2=list(V=1, nu=0.002)))

#Compute mean-scaled measurement error variances
mev=moddat$dse2/(moddat$mean^2)*100
moddat$SE=sqrt(mev)  

#Set sampling parameters
samples = 1000
thin = 10
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

#Sample MCMC
mod<-MCMCglmm(d ~ evals + tg1,
              random = ~study_ID + species,
              rcov = ~SE:units,
              mev = mev,
              data = moddat, 
              family = "gaussian", prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

#Check convergence
x11()
plot(mod$VCV)

summary(mod)

