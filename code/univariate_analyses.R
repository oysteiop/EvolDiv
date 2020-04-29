###############################################
#### - Univariate analysis of D database - ####
###############################################
rm(list=ls())

library(plyr)
library(reshape2)
library(lme4)
library(MCMCglmm)
library(devtools)
#install_github("GHBolstad/evolvability")
#library(withr)
#with_libpaths(new="C:/Program Files/R/R-3.5.0/library",install_github("GHBolstad/evolvability"))
#.libPaths("C:/Program Files/R/R-3.5.0/library")
library(evolvability)

ddat = read.table("data/dmatdata.txt", header=T)
ddat$ID = paste(ddat$reference, ddat$species, ddat$environment, sep="_")
maxdists = read.csv(file="data/maxdists.csv")

#Scaling of trait variances with the mean
#ddat=ddat[ddat$trait!="organ_size",]
#plot(log10(ddat$mean), log10(ddat$sd))
#lines(-10:10, -10:10)

#List of studies
studies = sort(unique(ddat$ID))
studies
splist = tapply(as.character(ddat$species), ddat$ID, function(x) x[1])

#Compute among-pop variance for each trait in each study
outlist = list()
for(s in 1:length(studies)){
  red = ddat[ddat$ID==studies[[s]],]
  red$trait = factor(red$trait)
  df=data.frame(study_ID = studies[[s]], 
                species = splist[[s]], 
                trait = sort(unique(red$trait)),
                environment = sort(unique(red$environment)),
                npop = tapply(red$mean>-Inf, red$trait, sum, na.rm=T),
                d = cbind(tapply(log(abs(red$mean)), red$trait, var, na.rm=T)),
                mean = cbind(tapply(red$mean, red$trait, mean, na.rm=T)))
  df = na.omit(df)
  outlist[[s]] = df
}
ddf = rbind.fill(outlist)
head(ddf, 5)

ddat$species_trait = paste0(ddat$species, "_", ddat$trait)
ddf$species_trait = paste0(ddf$species, "_", ddf$trait)

# Add metadata
tg1=NULL
for(i in 1:nrow(ddf)){
  tg1[i]=as.character(ddat$tg1)[which(as.character(ddat$species_trait)==as.character(ddf$species_trait)[i])[1]]
}

tg2=NULL
for(i in 1:nrow(ddf)){
  tg2[i]=as.character(ddat$tg2)[which(as.character(ddat$species_trait)==as.character(ddf$species_trait)[i])[1]]
}

dimension=NULL
for(i in 1:nrow(ddf)){
  dimension[i]=as.character(ddat$dimension)[which(as.character(ddat$species_trait)==as.character(ddf$species_trait)[i])[1]]
}

ms=NULL
for(i in 1:nrow(ddf)){
  ms[i]=as.character(ddat$ms)[which(as.character(ddat$species_trait)==as.character(ddf$species_trait)[i])[1]]
}

maxdist=NULL
for(i in 1:nrow(ddf)){
  maxdist[i]=as.numeric(maxdists[which(maxdists[,2]==as.character(ddf$study_ID)[i])[1],3])
}

ddf$maxdist = maxdist
ddf$tg1 = factor(tg1)
ddf$tg2 = factor(tg2)
ddf$ms = ms
ddf$ms = factor(ddf$ms, levels=c("S","M","O"))
ddf$dimension = dimension
ddf$dimension = factor(ddf$dimension, levels=c("linear", "area", "mass_volume", "count", "ratio", "time"))
ddf$environment = factor(ddf$environment, levels=c("greenhouse", "common_garden", "field"))

# Sampling variance of a variance (From Lynch & Walsh 1998 p. 815)
sv = (2*(ddf$d^2))/(ddf$npop+2)
ddf$d_se = sqrt(sv)
#plot(ddf$d, ddf$d_se)
#lines(0:100, 0:100)
#plot(ddf$npop, (ddf$d_se/ddf$d))

# Combine with data from Evolvability database
edat = read.table("data/evolvabilitydatabase2020.txt", header=T)
edat$species_measurement = paste0(edat$species,"_",edat$measurement)

evals = NULL
for(i in 1:nrow(ddf)){
w = which(as.character(edat$species)==as.character(ddf$species)[i] & as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])
evals[i] = mean(edat$evolvability[w], na.rm=T)
}
ddf$evals = evals

head(ddf)
summary(ddf)


red=ddf[ddf$ms%in%c("S","M","O"),]
red$ms=factor(red$ms,levels=c("S","M","O"))
red=red[red$tg1%in%c("floral","vegetative"),]
mstg=paste0(red$ms, red$tg1)
mstg=factor(mstg, levels=c("Sfloral","Svegetative","Mfloral","Mvegetative","Ofloral","Ovegetative"))
plot(factor(mstg), log(red$d),at=c(1,2,4,5,7,8), las=1, col=rep(c("white","grey"),2), xaxt="n", xlab="",
     ylab="Divergence", ylim=c(-15, 6))
legend("topleft", pch=c(15, 0), col=c("grey", "white"), legend=c("Vegetative", "Floral"), bty="n", pt.cex = 1.3)
legend("topleft", pch=c(0, 0), col=c("black", "black"), legend=c("", ""), bty="n", pt.cex=1.3)
axis(1, at=c(1.5, 4.5, 7.5), labels=c("Selfer", "Mixed", "Outcrossing"))

#### Summary stats with all divergence data####
length(unique(ddf$study_ID))
length(unique(ddf$species))

tapply(ddf$d*100, ddf$tg1, median)
round(sqrt(tapply(ddf$d, ddf$tg1, median))*sqrt(2/pi)*100, 1) #dL
tapply(ddf$d>-Inf, ddf$tg1, sum)

tapply(ddf$d*100, list(ddf$tg1, ddf$tg2), median)
round(sqrt(tapply(ddf$d, list(ddf$tg1, ddf$tg2), median))*sqrt(2/pi)*100, 1) #dL
tapply(ddf$d>-Inf, ddf$tg2, sum)

tapply(ddf$d*100, list(ddf$tg1, ddf$dimension), median)

tapply(ddf$d*100, list(ddf$ms, ddf$tg1), median, na.rm=T)
tapply(ddf$d>-Inf, list(ddf$ms, ddf$tg1), sum, na.rm=T)

lin = ddf[ddf$dimension=="linear",]
tapply(lin$d*100, lin$tg1, median)
t(tapply(lin$d*100, list(lin$tg1,lin$species), median))

are = ddf[ddf$dimension=="area",]
t(tapply(are$d*100, list(are$tg1,are$species), median))

mvo = ddf[ddf$dimension=="mass_volume",]
t(tapply(mvo$d*100, list(mvo$tg1,mvo$species), median))

cou = ddf[ddf$dimension=="count",]
t(tapply(cou$d*100, list(cou$tg1, cou$species), median))

# Dimensions plot
medians = tapply(log(ddf$evals), ddf$dimension, median, na.rm=T)
ddf$dimension = factor(ddf$dimension, levels=names(sort(medians, decreasing=F)))

x11(height=4.5, width=4.5)
plot(ddf$dimension, log10(ddf$d*100), at=c(1,4,7,10,13,16), xlim=c(0,17), ylim=c(-3, 5.5), col="grey", xaxt="n", las=1, yaxt="n")
par(new=T)
plot(ddf$dimension, log10(ddf$evals), at=c(2,5,8,11,14,17), xlim=c(0,17), ylim=c(-3, 5.5), xaxt="n", yaxt="n")
axis(1, at=c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5), labels=rep("", 6))
text(c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5), par("usr")[3] - .5, srt = 45, adj = 1,cex=1,
     labels = levels(ddf$dimension), xpd = TRUE)
legend("topleft", pch=c(15, 0), col=c("grey", "white"), legend=c("Evolvability (%)", "Divergence (x100)"), bty="n", pt.cex = 1.2)
legend("topleft", pch=c(22, 0), col=c("black", "black"), legend=c("", ""), bty="n", pt.cex=1.4)
axis(2, at=c(-2:3), signif(10^(-2:3),1), las=1)


# Boxplots ####
ddf$logd = log10(ddf$d*100)
ddf$logd[which(ddf$logd==-Inf)]=log10(0.001) #Set zero evolvabilities to 0.01
ddf$logd[which(ddf$logd<(-3))]=log10(0.001) #Set zero evolvabilities to 0.01

tpos=-.25

veg = ddf[ddf$tg1=="vegetative",]

medians = tapply(veg$logd, veg$tg2, median, na.rm=T)
veg$tg2 = factor(veg$tg2, levels=names(sort(medians, decreasing=T)))
levels(veg$tg2)
a = 1:7

x11(height=5, width=6.5)
par(mar=c(6,4,2,2))
plot(veg$tg2, veg$logd, cex=.8, notch=F, xaxt="n", col="darkgreen", 
     ylab="Divergence (x100)", xlab="", yaxt="n", ylim=c(-3,3), xlim=c(0,30))
axis(2, at=c(-3:3), c("<0.001",c(signif(10^(c(-2:3)),1))), las=1)
axis(1, at=a, labels = FALSE)
labels = paste0(toupper(as.character(levels(veg$tg2)))," (",tapply(veg$d>-100, veg$tg2, sum, na.rm=T),") ")
labels = sub("_"," ",labels)
text(a, par("usr")[3] + tpos, srt = 45, adj = 1,cex=.7,
     labels = labels, xpd = TRUE)
h = median(ddf$logd[ddf$tg1=="vegetative"], na.rm=T)
segments(min(a), h, max(a), h, lwd=3)

#####

lifehist = ddf[ddf$tg1=="lifehistory",]
medians = tapply(lifehist$logd, lifehist$tg2, median, na.rm=T)
lifehist$tg2 = factor(lifehist$tg2, levels=names(sort(medians,decreasing=T)))
levels(lifehist$tg2)
a = (max(a)+2):((max(a)+1)+length(levels(lifehist$tg2)))
par(new=T)
plot(lifehist$tg2, lifehist$logd, cex=.8, yaxt="n", at=a, xaxt="n", col="brown",
     ylab="", xlab="", ylim=c(-3,3), xlim=c(0,30))
axis(1, at=a, labels = FALSE)
labels = paste0(toupper(as.character(levels(lifehist$tg2)))," (",tapply(lifehist$logd>-100, lifehist$tg2, sum, na.rm=T),") ")
labels = sub("_"," ",labels)
text(a, par("usr")[3] + tpos, srt = 45, adj = 1, cex=.7,
     labels = labels, xpd = TRUE)
h = median(ddf$logd[ddf$tg1=="lifehistory"], na.rm=T)
segments(min(a), h, max(a), h, lwd=3)

####

floral = ddf[ddf$tg1=="floral" & ddf$tg2!="fitness",]
#floral$logd[which(floral$logd==-Inf)]=log10(0.01)

medians = tapply(floral$logd, floral$tg2, median, na.rm=T)
floral$tg2 = factor(floral$tg2, levels=names(sort(medians,decreasing=T)))
levels(floral$tg2)
a = (max(a)+2):((max(a)+1)+length(levels(floral$tg2)))

par(new=T)
plot(floral$tg2, floral$logd, cex=.8, yaxt="n", at=a, xaxt="n", col="darkblue",
     ylab="", xlab="", ylim=c(-3,3), xlim=c(0,30))
axis(1, at=a, labels = FALSE)
labels = paste0(toupper(as.character(levels(floral$tg2)))," (",tapply(floral$logd>-100, floral$tg2, sum, na.rm=T),") ")
labels = sub("_"," ",labels)
text(a, par("usr")[3] + tpos, srt = 45, adj = 1, cex=.7,
     labels = labels, xpd = TRUE)
h = median(ddf$logd[ddf$tg1=="floral"], na.rm=T)
segments(min(a), h, max(a), h, lwd=3)

#### Subset data for d vs. e meta-analysis ####
ddf = na.omit(ddf)
ddf = ddf[ddf$d>0,]
ddf = ddf[ddf$evals>0,]
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

#### Summary stats ####
length(unique(ddf$study_ID))
length(unique(ddf$species))

tapply(ddf$d*100, ddf$tg1, median)
tapply(ddf$d>-Inf, ddf$tg1, sum)
tapply(ddf$d*100, list(ddf$tg1, ddf$tg2), median)
tapply(ddf$d>-Inf, ddf$tg2, sum)

tapply(ddf$d*100, list(ddf$tg1, ddf$dimension), median)

tapply(ddf$d*100, list(ddf$ms, ddf$tg1), median, na.rm=T)
tapply(ddf$d>-Inf, list(ddf$ms, ddf$tg1), sum, na.rm=T)

lin = ddf[ddf$dimension=="linear",]
t(tapply(lin$d*100, list(lin$tg1,lin$species), median))

#### Informal meta-analysis ####
mdat = na.omit(subset(lin[lin$d>0,], select=c("d","evals","npop","maxdist","species","study_ID")))
names(mdat)
m0 = lmer(log(d)~ log(evals) + log(npop) + log(maxdist) + (1|species/study_ID), REML=T, data=mdat)
summary(m0)

mdat = mdat[mdat$study_ID!="Hansen_et_al._2003_Dalechampia_scandens_A_greenhouse",]
mdat = mdat[mdat$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_field",]
m0 = lmer(log(d)~ log(evals) + log(npop) + log(maxdist) + (1|species/study_ID), REML=T, data=mdat)
summary(m0)


# Remove some of the repeated D. scandens studies?
ddf = ddf[ddf$study_ID!="Hansen_et_al._2003_Dalechampia_scandens_A_greenhouse",]
ddf = ddf[ddf$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_field",]
#ddf=ddf[ddf$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_greenhouse",]

ddf$scale_npop = scale(log(ddf$npop), scale=F)
ddf$scale_maxdist = scale(log(ddf$maxdist), scale=F)
wgts = 1/((ddf$d_se^2)/(ddf$d^2))
plot(ddf$npop, wgts)

m0 = lmer(log(d)~ log(evals) + log(npop) + log(maxdist) + (log(evals)|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
summary(m0)

m0b = lmer(log(d)~ log(evals) + log(npop) + log(maxdist) + (log(evals)|species/study_ID), REML=F, data=ddf)
AIC(m0,m0b)
summary(m0b)

m = lmer(log(d)~ log(evals)*log(maxdist) + log(npop) + (log(evals)|study_ID) +(1|species), REML=F, weights=wgts, data=ddf)
summary(m)
AIC(m0, m)

m = lmer(log(d)~ log(evals) + log(npop) + log(maxdist) + (log(evals)|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
summary(m)

m = lmer(log(d)~ log(evals)*dimension + log(npop) + log(maxdist) + (log(evals)|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
logLik(m0)
logLik(m)

summary(m)

m = lmer(log(d)~ log(evals)*environment + log(npop) + log(maxdist) + (log(evals)|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
summary(m)

floveg = ddf[ddf$tg1=="floral" | ddf$tg1=="vegetative",]
floveg$tg1 = factor(floveg$tg1)
m0 = lmer(log(d) ~ log(evals)+ log(npop) +log(maxdist) + (1|species/study_ID), REML=F, data=floveg)
m = lmer(log(d) ~ log(evals)*tg1 + log(npop) + log(maxdist) + (1|species/study_ID), REML=F, data=floveg)
m = lmer(log(d) ~ log(evals)*log(maxdist) + log(npop) + (1|species/study_ID), REML=F, data=floveg)

AIC(m0, m)
AIC(m0, m)[2,2]-AIC(m0, m)[1,2]

m = lmer(log(d) ~ -1+ tg1 + log(evals):tg1 + log(npop) + log(maxdist) + (1|species/study_ID), data=floveg)
summary(m)





# Formal meta-analysis using Almer_SE
SE = sqrt((ddf$d_se^2)/(ddf$d^2))
plot(ddf$npop, SE)
m0 = evolvability::Almer_SE(log(d) ~ log(evals) + log(npop) + log(maxdist) + (log(evals)|study_ID), 
           SE=SE, maxiter = 100, data=ddf)
m1 = evolvability::Almer_SE(log(d) ~ log(evals) + log(npop) + log(maxdist) + (log(evals)|study_ID) + (1|species), 
                            SE=SE, maxiter = 100, data=ddf)

AIC(m0, m1)
summary(m1)

#### Plotting evolvability vs. divergence ####

#ddf$d[which(ddf$d<0.00001)]=0.00001

#All data
par(mfrow=c(1,1))
plot(log10(ddf$evals),ddf$logd,
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf$npop),
     col="black",
     xlim=c(-2.5, 2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

#Floral and vegetative
floveg = ddf[ddf$tg1=="floral" | ddf$tg1=="vegetative",]
floveg$tg1 = factor(floveg$tg1)
m = lmer(log(d*100) ~ -1 + tg1 + log(evals):tg1 + scale_npop+ scale_maxdist 
         + (1|species/study_ID), data=floveg)

x11(height=5, width=5)
par(mfrow=c(1,1))
plot(log10(floveg$evals), floveg$logd,
     xlab="Evolvability (%)",
     ylab="Divergence (x100)",
     pch=1, cex=1*sqrt(ddf$npop),
     col=c(rgb(0, 0, 0.545, .35), 
           rgb(0.004, 0.196, 0.125, .35))[as.numeric(floveg$tg1)],
     xlim=c(-2.5, 2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

x1=seq(min(log(floveg$evals[floveg$tg1=="floral"])),
       max(log(floveg$evals[floveg$tg1=="floral"])), .1)
y = summary(m)$coef[1,1] + summary(m)$coef[5,1]*x1 + 
    summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="floral"]) +
    summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="floral"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e = log(median(ddf$evals[ddf$tg1=="floral"]))
y = summary(m)$coef[1,1] + summary(m)$coef[5,1]*me_e +
  summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="floral"]) +
  summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="floral"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(floveg$evals[floveg$tg1=="vegetative"])),
       max(log(floveg$evals[floveg$tg1=="vegetative"])), .1)
y = summary(m)$coef[2,1] + summary(m)$coef[6,1]*x2 +
    summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="vegetative"]) +
    summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="vegetative"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$tg1=="vegetative"]))
y = summary(m)$coef[2,1] + summary(m)$coef[6,1]*me_e +
  summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="vegetative"]) +
  summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="vegetative"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

legend("topleft", pch=15, col=c("darkblue", "darkgreen"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Floral", "Vegetative"))

####

x11(height=5, width=6)
par(mfrow=c(1,1))
plot(log10(floveg$evals), floveg$logd,
     xlab="Evolvability (%)",
     ylab="Divergence (x100)",
     pch=1, cex=1*sqrt(ddf$npop),
     col=c(rgb(0, 0, 0.545, .35), 
           rgb(0.004, 0.196, 0.125, .35))[as.numeric(floveg$tg1)],
     xlim=c(-2.5, 3.5), ylim=c(-3.5,3), xaxt="n", yaxt="n")
axis(1, c(-2,-1,0,1,2), 10^c(-2,-1,0,1,2))
axis(2, c(-3,-2,-1,0,1,2), c("<0.001", 10^c(-2,-1,0,1,2)), las=1)

par(new=T)
plot(floveg$tg1, floveg$logd, at=c(2.5, 3), boxwex=0.4, xlab="", ylab="",
     col=c(rgb(0, 0, 0.545, .35), rgb(0.004, 0.196, 0.125, .35)),
     xlim=c(-2.5, 3.5), ylim=c(-3.5,3), xaxt="n", yaxt="n")

x1=seq(min(log(floveg$evals[floveg$tg1=="floral"])),
       max(log(floveg$evals[floveg$tg1=="floral"])), .1)
y = summary(m)$coef[1,1] + summary(m)$coef[5,1]*x1 + 
  summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="floral"]) +
  summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="floral"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e = log(median(ddf$evals[ddf$tg1=="floral"]))
y = summary(m)$coef[1,1] + summary(m)$coef[5,1]*me_e +
  summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="floral"]) +
  summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="floral"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(floveg$evals[floveg$tg1=="vegetative"])),
       max(log(floveg$evals[floveg$tg1=="vegetative"])), .1)
y = summary(m)$coef[2,1] + summary(m)$coef[6,1]*x2 +
  summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="vegetative"]) +
  summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="vegetative"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$tg1=="vegetative"]))
y = summary(m)$coef[2,1] + summary(m)$coef[6,1]*me_e +
  summary(m)$coef[3,1]*mean(ddf$scale_npop[ddf$tg1=="vegetative"]) +
  summary(m)$coef[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="vegetative"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

legend("topleft", pch=15, col=c("darkblue", "darkgreen"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Floral", "Vegetative"))


# Mating systems
m = lmer(log(d*100) ~ -1 + ms + log(evals):ms + scale_npop+ scale_maxdist  
         + (1|species/study_ID), data=ddf)
summary(m)$coef

x11(height=5, width=6)
par(mfrow=c(1,1))
plot(log10(ddf$evals), ddf$logd,
     xlab="",
     ylab="Divergence (x100)",
     pch=1, cex=1*sqrt(ddf$npop),
     col=c(rgb(0, 0, 0.545, .35), 
           rgb(0.004, 0.196, 0.125, .35),
           rgb(.545, 0, 0, .35))[as.numeric(ddf$ms)],
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")
axis(1, c(-2,-1,0,1,2), 10^c(-2,-1,0,1,2))
axis(2, c(-3,-2,-1,0,1,2), c("<0.001", 10^c(-2,-1,0,1,2)), las=1)
mtext("Evolvability (%)", 1, line=2.5)

par(new=T)
plot(ddf$ms, ddf$logd, at=c(2.5, 3, 3.5), boxwex=0.4, xlab="", ylab="",
     col=c(rgb(0, 0, 0.545, .35), rgb(0.004, 0.196, 0.125, .35), rgb(.545, 0, 0, .35)),
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")

x1=seq(min(log(ddf$evals[ddf$ms=="S"])),
       max(log(ddf$evals[ddf$ms=="S"])), .1)
y = summary(m)$coef[1,1] + summary(m)$coef[6,1]*x1 +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$ms=="S"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$ms=="S"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e=log(median(ddf$evals[ddf$ms=="S"]))
y = summary(m)$coef[1,1] + summary(m)$coef[6,1]*me_e +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$ms=="S"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$ms=="S"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(ddf$evals[ddf$ms=="M"])),
       max(log(ddf$evals[ddf$ms=="M"])), .1)
y = summary(m)$coef[2,1] + summary(m)$coef[7,1]*x2 +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$ms=="M"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$ms=="M"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$ms=="M"])) 
y = summary(m)$coef[2,1] + summary(m)$coef[7,1]*me_e +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$ms=="M"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$ms=="M"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

x3=seq(min(log(ddf$evals[ddf$ms=="O"])),
       max(log(ddf$evals[ddf$ms=="O"])), .1)
y = summary(m)$coef[3,1] + summary(m)$coef[8,1]*x3 +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$ms=="O"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$ms=="O"])
lines(log10(exp(x3)), log10(exp(y)), col ="darkred", lwd=2)

me_e = log(median(ddf$evals[ddf$ms=="O"]))
y = summary(m)$coef[3,1] + summary(m)$coef[8,1]*me_e +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$ms=="O"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$ms=="O"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkred")

legend("topleft", pch=15, col=c("darkblue", "darkgreen", "darkred"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Selfing", "Mixed", "Outcrossing"))

# Environments
m = lmer(log(d*100) ~ -1 + environment + log(evals):environment + scale_npop+ scale_maxdist  
         + (1|species/study_ID), data=ddf)
summary(m)$coef

x11(height=5, width=6)
par(mfrow=c(1,1))
plot(log10(ddf$evals), ddf$logd,
     xlab="",
     ylab="Divergence (x100)",
     pch=1, cex=1*sqrt(ddf$npop),
     col=c(rgb(0, 0, 0.545, .35), 
           rgb(0.004, 0.196, 0.125, .35),
           rgb(.545, 0, 0, .35))[as.numeric(ddf$environment)],
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")
axis(1, c(-2,-1,0,1,2), 10^c(-2,-1,0,1,2))
axis(2, c(-3,-2,-1,0,1,2), c("<0.001", 10^c(-2,-1,0,1,2)), las=1)
mtext("Evolvability (%)", 1, line=2.5)

par(new=T)
plot(ddf$environment, ddf$logd, at=c(2.5, 3, 3.5), boxwex=0.4, xlab="", ylab="",
     col=c(rgb(0, 0, 0.545, .35), rgb(0.004, 0.196, 0.125, .35), rgb(.545, 0, 0, .35)),
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")

x1=seq(min(log(ddf$evals[ddf$environment=="greenhouse"])),
       max(log(ddf$evals[ddf$environment=="greenhouse"])), .1)
y = summary(m)$coef[1,1] + summary(m)$coef[6,1]*x1 +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$environment=="greenhouse"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$environment=="greenhouse"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e=log(median(ddf$evals[ddf$environment=="greenhouse"]))
y = summary(m)$coef[1,1] + summary(m)$coef[6,1]*me_e +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$environment=="greenhouse"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$environment=="greenhouse"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(ddf$evals[ddf$environment=="common_garden"])),
       max(log(ddf$evals[ddf$environment=="common_garden"])), .1)
y = summary(m)$coef[2,1] + summary(m)$coef[7,1]*x2 +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$environment=="common_garden"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$environment=="common_garden"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$environment=="common_garden"])) 
y = summary(m)$coef[2,1] + summary(m)$coef[7,1]*me_e +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$environment=="common_garden"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$environment=="common_garden"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

x3=seq(min(log(ddf$evals[ddf$environment=="field"])),
       max(log(ddf$evals[ddf$environment=="field"])), .1)
y = summary(m)$coef[3,1] + summary(m)$coef[8,1]*x3 +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$environment=="field"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$environment=="field"])
lines(log10(exp(x3)), log10(exp(y)), col ="darkred", lwd=2)

me_e = log(median(ddf$evals[ddf$environment=="field"]))
y = summary(m)$coef[3,1] + summary(m)$coef[8,1]*me_e +
  summary(m)$coef[4,1]*mean(ddf$scale_npop[ddf$environment=="field"]) +
  summary(m)$coef[5,1]*mean(ddf$scale_maxdist[ddf$environment=="field"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkred")

legend("topleft", pch=15, col=c("darkblue", "darkgreen", "darkred"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Greenhouse", "Garden", "Field"))


# Dimensions
library(fBasics)
m = lmer(log(d*100) ~ -1 + dimension + log(evals):dimension + scale_npop+ scale_maxdist 
         + (1|species/study_ID), data=ddf)
cols = c("blue3", "firebrick", "green3", "yellow3","cyan", "black") 

x11(height=4.5, width=6.5)
par(mfrow=c(1,1), mar=c(4,4,2,7))
plot(log10(ddf$evals), log10(ddf$d*100),
     xlab="Evolvability (%)",
     ylab="Divergence (x100)",
     pch=1, cex=1*sqrt(ddf$npop),
     col=cols[as.numeric(ddf$dimension)],
     xlim=c(-2.5, 2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

x1=seq(min(log(ddf$evals[ddf$dimension=="ratio"])),
       max(log(ddf$evals[ddf$dimension=="ratio"])), .1)
y = summary(m)$coef[1,1] + summary(m)$coef[9,1]*x1 + 
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="ratio"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="ratio"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[1], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="linear"])),
       max(log(ddf$evals[ddf$dimension=="linear"])), .1)
y = summary(m)$coef[2,1] + summary(m)$coef[10,1]*x1 + 
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="linear"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="linear"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[2], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="time"])),
       max(log(ddf$evals[ddf$dimension=="time"])), .1)
y = summary(m)$coef[3,1] + summary(m)$coef[11,1]*x1 + 
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="time"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="time"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[3], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="mass_volume"])),
       max(log(ddf$evals[ddf$dimension=="mass_volume"])), .1)
y = summary(m)$coef[4,1] + summary(m)$coef[12,1]*x1 + 
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="mass_volume"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="mass_volume"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[4], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="area"])),
       max(log(ddf$evals[ddf$dimension=="area"])), .1)
y = summary(m)$coef[5,1] + summary(m)$coef[13,1]*x1 + 
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="area"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="area"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[5], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="count"])),
       max(log(ddf$evals[ddf$dimension=="count"])), .1)
y = summary(m)$coef[6,1] + summary(m)$coef[14,1]*x1 + 
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="count"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="count"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[6], lwd=2)

me_e = log(median(ddf$evals[ddf$dimension=="ratio"]))
y = summary(m)$coef[1,1] + summary(m)$coef[9,1]*me_e +
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="ratio"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="ratio"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="linear"]))
y = summary(m)$coef[2,1] + summary(m)$coef[10,1]*me_e +
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="linear"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="linear"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="time"]))
y = summary(m)$coef[3,1] + summary(m)$coef[11,1]*me_e +
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="time"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="time"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="mass_volume"]))
y = summary(m)$coef[4,1] + summary(m)$coef[12,1]*me_e +
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="mass_volume"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="mass_volume"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="area"]))
y = summary(m)$coef[5,1] + summary(m)$coef[13,1]*me_e +
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="area"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="area"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="count"]))
y = summary(m)$coef[6,1] + summary(m)$coef[14,1]*me_e +
  summary(m)$coef[7,1]*mean(ddf$scale_npop[ddf$dimension=="count"]) +
  summary(m)$coef[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="count"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

legend(x=2.2, y = 2, pch=15, col=cols, cex=1, pt.cex=1.5, 
       bty="n", legend=levels(ddf$dimension), xpd=T)


# Contextual model
dimmeans = tapply(log(ddf$evals), ddf$dimension, mean)

meanvals=NULL
for(i in 1:nrow(ddf)){
 meanvals[i]=dimmeans[which(names(dimmeans)==ddf$dimension[i])] 
}
data.frame(meanvals, ddf$dimension)

m = lmer(log(d*100) ~ log(evals) + meanvals + scale_npop + scale_maxdist
         + (1|species/study_ID), data=ddf)
summary(m)


# Removing the dimension effect
resd = summary(lm(log10(ddf$d*100)~ddf$dimension))$residuals
rese = summary(lm(log10(ddf$evals)~ddf$dimension))$residuals

plot(rese, resd)
m = lmer(resd ~ rese + scale_npop + scale_maxdist
         + (1|species/study_ID), data=ddf)
summary(m)

x11(height = 4.5, width=6.5)
par(mfrow=c(1,1), mar=c(4,4,2,7))
plot(rese, resd, las=1,
     xlab="Residual evolvability (%)",
     ylab="Residual divergence (x100)",
     pch=1, cex=1*sqrt(ddf$npop),
     col=cols[as.numeric(ddf$dimension)])

#axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
#axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)
legend(x=2.4, y = 2, pch=15, col=cols, cex=1, pt.cex=1.5, 
       bty="n", legend=levels(ddf$dimension), xpd=T)

# Function to plot e vs. d with subset highlighted ####
plotSubset = function(category, subset, xlim=c(-2.5,2), ylim=c(-4,3), ...){
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
par(mfrow=c(1,3))
plotSubset("environment", "greenhouse", lwd=2)
plotSubset("environment", "field", lwd=2)
plotSubset("environment", "common_garden", lwd=2)

x11()
par(mfrow=c(1,2))
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
plotSubset("dimension", "time")


studies = unique(ddf$study_ID)

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
mev = moddat$d_se^2
moddat$SE = sqrt(mev)  

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
str(summary(mod)$solutions)
ests = summary(mod)$solutions

x11(height=5, width=5)
par(mar=c(4, 8, 2, 2))
plot(rev(ests[-1,1]), 1:12, xlim=c(-4, 3), las=1, pch=16, yaxt="n", ylab="", xlab="")
mtext(1, text="Estimate (95% CI)", line=2.5)
axis(2, 1:12, labels=F)
segments(rev(ests[-1,2]), 1:12, rev(ests[-1,3]), 1:12)
abline(v=0, lty=2)
labels=rev(rownames(ests)[-1])
text(x=-4.7, y=1:12, labels=labels, las=1, cex=.8, adj=1, xpd=T)


plot(log(moddat$maxdist), log(moddat$d))
plot(log(moddat$evals), log(moddat$d))
plot(log(moddat$npop), log(moddat$d))


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
