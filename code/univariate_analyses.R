###############################################
#### - Univariate analysis of D database - ####
###############################################
rm(list=ls())

library(devtools)
#install_github("GHBolstad/evolvability")
#library(withr)
#with_libpaths(new="C:/Program Files/R/R-3.5.0/library",install_github("GHBolstad/evolvability"))
#.libPaths("C:/Program Files/R/R-3.5.0/library")
library(plyr)
library(reshape2)
library(lme4)
library(glmmTMB)
library(MCMCglmm)
library(evolvability)
library(MuMIn)
library(fBasics)

nbot = 10

# Read datafiles
ddat = read.table("data/dmatdata.txt", header=T)
ddat$ID = paste(ddat$reference, ddat$species, ddat$environment, sep="_")
maxdists = read.csv(file="data/maxdists.csv")

# List of studies
studies = sort(unique(ddat$ID))
studies
splist = tapply(as.character(ddat$species), ddat$ID, function(x) x[1])

# Compute among-pop variance for each trait in each study
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
                mean = cbind(tapply(red$mean, red$trait, mean, na.rm=T)),
                ref = unique(red$reference))
  df = na.omit(df)
  outlist[[s]] = df
}
ddf = rbind.fill(outlist)
head(ddf, 5)

ddat$species_trait = paste0(ddat$species, "_", ddat$trait)
ddf$species_trait = paste0(ddf$species, "_", ddf$trait)

length(unique(ddf$species))
length(unique(ddf$species_trait))
length(unique(ddat$population))

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

tmvals = NULL
for(i in 1:nrow(ddf)){
  w = which(as.character(edat$species)==as.character(ddf$species)[i] & as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])
  tmvals[i] = mean(edat$outcrossing_best[w], na.rm=T)
}
#ddf$tmvals = tmvals

# Study table
stab = ddply(ddf, .(study_ID), summarize,
             species = unique(species),
             npop = median(npop),
             ntraits = length(unique(trait)),
             environment = unique(environment),
             ref = unique(ref))
head(stab)
write.csv2(stab[order(stab$species),], "studytable.csv")


#### Start of analyses ####
head(ddf)
summary(ddf)

#### Summary stats with all divergence data ####
length(unique(ddf$study_ID)) # Number of divergence studies
length(!is.na(ddat$mean)) # Number of population means
length(unique(paste0(ddat$species, ddat$population))) # Number of populations
length(unique(ddf$species)) # Number of species
length(unique(ddat$family)) # Number of families

ddf$dP = exp(sqrt(ddf$d)*sqrt(2/pi)) # Compute proportional divergence metric

tapply(ddf$d, ddf$tg1, median)
tapply(ddf$dP, ddf$tg1, median)
tapply(ddf$d>-Inf, ddf$tg1, sum)

# Standard errors for floral and vegetative divergence
flosamp = NULL
for(i in 1:nbot){
flosamp[i] = median(sample(ddf$dP[which(ddf$tg1=="floral")], replace=T))
}
vegsamp = NULL
for(i in 1:nbot){
  vegsamp[i] = median(sample(ddf$dP[which(ddf$tg1=="vegetative")], replace=T))
}

sd(vegsamp)
sd(flosamp)

# Model comparison flo vs. veg
floveg = ddf[ddf$tg1=="floral" | ddf$tg1=="vegetative",]
floveg$tg1 = factor(floveg$tg1)
floveg = floveg[floveg$d>0,]

m0 = glmmTMB(log(d)~ 1 + (1|study_ID) + (1|species), REML=F, data=floveg)
m = glmmTMB(log(d)~ tg1 + (1|study_ID) + (1|species), REML=F, data=floveg)

AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1] #delta AIC, support for model m

summary(m)

#floveg = floveg[!is.na(floveg$ms),]

#m0 = lmer(log(d)~ 1 + (1|study_ID) + (1|species), REML=F, data=floveg)
#m1 = lmer(log(d)~ tg1 + (1|study_ID) + (1|species), REML=F, data=floveg)
#m2 = lmer(log(d)~ ms + (1|study_ID) + (1|species), REML=F, data=floveg)
#m3 = lmer(log(d)~ tg1 + ms + (1|study_ID) + (1|species), REML=F, data=floveg)
#m4 = lmer(log(d)~ tg1 * ms + (1|study_ID) + (1|species), REML=F, data=floveg)

#AIC(m0, m1, m2, m3, m4)

# Linear traits only
lin = ddf[ddf$dimension=="linear",]
tapply(lin$d, lin$tg1, median)
tapply(lin$dP, lin$tg1, median)

flosamp = NULL
for(i in 1:nbot){
  flosamp[i] = median(sample(lin$dP[which(lin$tg1=="floral")], replace=T))
}
vegsamp = NULL
for(i in 1:nbot){
  vegsamp[i] = median(sample(lin$dP[which(lin$tg1=="vegetative")], replace=T))
}

sd(vegsamp)
sd(flosamp)

# Linear floral vs. linear vegetative
floveglin = floveg[floveg$dimension=="linear",]

m0 = glmmTMB(log(d)~ 1 + (1|study_ID) + (1|species), REML=F, data=floveglin)
m = glmmTMB(log(d)~ tg1 + (1|study_ID) + (1|species), REML=F, data=floveglin)

AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1]

summary(m)

# By study
pairwise = na.omit(as.data.frame(tapply(lin$dP, list(lin$study_ID, lin$tg1), median, na.rm=T)[,c(1,4)]))
pairwise
sum(pairwise$vegetative>pairwise$floral)
sum(!is.na(pairwise$vegetative>pairwise$floral))

# Subcategories Within floral traits
flo = droplevels(ddf[ddf$tg1=="floral",])
sort(tapply(flo$d, flo$tg2, median, na.rm=T))
sort(tapply(flo$dP, flo$tg2, median, na.rm=T))
sort(tapply(flo$dP>-1, flo$tg2, sum, na.rm=T))

# Summary stats per mating system category ####
tapply(ddf$d, ddf$ms, median)
tapply(ddf$dP, ddf$ms, median)
tapply(ddf$d>-Inf, ddf$ms, sum)

selfsamp = NULL
for(i in 1:nbot){
  selfsamp[i] = median(sample(ddf$dP[which(ddf$ms=="S")], replace=T))
}
mixsamp = NULL
for(i in 1:nbot){
  mixsamp[i] = median(sample(ddf$dP[which(ddf$ms=="M")], replace=T))
}
outsamp = NULL
for(i in 1:nbot){
  outsamp[i] = median(sample(ddf$dP[which(ddf$ms=="O")], replace=T))
}

sd(selfsamp)
sd(mixsamp)
sd(outsamp)
quantile(selfsamp, c(0.025, 0.975))
quantile(mixsamp, c(0.025, 0.975))
quantile(outsamp, c(0.025, 0.975))

#ss = na.omit(subset(ddf, select=c("d","d_se","ms")))
#ss = ss[ss$d>0,]
#ss$se2 = (ss$d_se^2)/(ss$d^2)
#summary(lm(log(d)~ms, weights = 1/se2, data=ss))
#plot(ss$ms, log(ss$d))

round(tapply(ddf$d, list(ddf$ms, ddf$tg1), median), 3)
round(tapply(ddf$d>0, list(ddf$ms, ddf$tg1), sum, na.rm=T), 3)

# Mating system vs. trait group ####
red = ddf[ddf$ms%in%c("S","M","O"),]
red$ms = factor(red$ms, levels=c("S","M","O"))
red = red[red$tg1%in%c("floral","vegetative"),]
mstg = paste0(red$ms, red$tg1)
mstg = factor(mstg, levels=c("Sfloral","Svegetative","Mfloral","Mvegetative","Ofloral","Ovegetative"))

#x11(height=4.5, width=4.5)
plot(factor(mstg), log10(red$d*100), at=c(1,2,4,5,7,8), las=1, col=rep(c("white","grey"),2), xaxt="n", xlab="",
     ylab="Proportional divergence", ylim=c(-3, 4), yaxt="n")
legend("topleft", pch=c(15, 0), col=c("grey", "white"), legend=c("Vegetative", "Floral"), bty="n", pt.cex = 1.3)
legend("topleft", pch=c(0, 0), col=c("black", "black"), legend=c("", ""), bty="n", pt.cex=1.3)
axis(1, at=c(1.5, 4.5, 7.5), labels=c("Selfer", "Mixed", "Outcrossing"))

xt3 = c(1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
x3at = log10(100*log(xt3)^2/(2/pi))
axis(2, at=x3at, signif(xt3, 4), las=1)


# Summary stats per study environment category ####
tapply(ddf$d, ddf$environment, median)
tapply(ddf$dP, ddf$environment, median)
tapply(ddf$d>-Inf, ddf$environment, sum)

ghsamp = NULL
for(i in 1:nbot){
  ghsamp[i] = median(sample(ddf$dP[which(ddf$environment=="greenhouse")], replace=T))
}
cgsamp = NULL
for(i in 1:nbot){
  cgsamp[i] = median(sample(ddf$dP[which(ddf$environment=="common_garden")], replace=T))
}
fsamp = NULL
for(i in 1:nbot){
  fsamp[i] = median(sample(ddf$dP[which(ddf$environment=="field")], replace=T))
}
sd(ghsamp)
sd(cgsamp)
sd(fsamp)

#ss = na.omit(subset(ddf, select=c("d","environment")))
#ss = ss[ss$d>0,]
#summary(lm(log(d)~environment, data=ss))
#plot(ss$environment, log(ss$d))

# Environment vs. trait group ####
tapply(ddf$d, list(ddf$tg1, ddf$environment), median)
tapply(ddf$dP, list(ddf$tg1, ddf$environment), median)
tapply(ddf$d>-Inf, list(ddf$tg1, ddf$environment), sum)

red = ddf[ddf$environment%in%c("field","greenhouse"),]
red$environment = factor(red$environment, levels=c("field","greenhouse"))
red = red[red$tg1%in%c("floral","vegetative"),]
mstg = paste0(red$environment, red$tg1)
mstg = factor(mstg, levels=c("fieldfloral","fieldvegetative","greenhousefloral","greenhousevegetative"))

#x11(height=4.5, width=4.5)
plot(factor(mstg), log10(red$d*100), at=c(1,2,4,5), las=1, col=rep(c("white","grey"),2), xaxt="n", xlab="",
     ylab="Proportional divergence", ylim=c(-3, 4), yaxt="n")
legend("topleft", pch=c(15, 0), col=c("grey", "white"), legend=c("Vegetative", "Floral"), bty="n", pt.cex = 1.3)
legend("topleft", pch=c(0, 0), col=c("black", "black"), legend=c("", ""), bty="n", pt.cex=1.3)
axis(1, at=c(1.5, 4.5), labels=c("Field", "Greenhouse"))

xt3 = c(1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
x3at = log10(100*log(xt3)^2/(2/pi))
axis(2, at=x3at, signif(xt3, 4), las=1)

# Trait dimensions plot ####
medians = tapply(log(ddf$evals), ddf$dimension, median, na.rm=T)
ddf$dimension = factor(ddf$dimension, levels=names(sort(medians, decreasing=F)))

#x11(height=4.5, width=4.5)
par(mar=c(5,4,2,4))
plot(ddf$dimension, log10(ddf$d*100), at=c(1,4,7,10,13,16), xlim=c(0,17), ylim=c(-3, 5.5), col="grey", xaxt="n", las=1, yaxt="n", xlab="")
par(new=T)
plot(ddf$dimension, log10(ddf$evals), at=c(2,5,8,11,14,17), xlim=c(0,17), ylim=c(-3, 5.5), xaxt="n", yaxt="n", xlab="")
axis(1, at=c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5), labels=rep("", 6))
text(c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5), par("usr")[3] - .5, srt = 45, adj = 1,cex=1,
     labels = levels(ddf$dimension), xpd = TRUE)
legend("topleft", pch=c(15, 0), col=c("grey", "white"), legend=c("Evolvability (%)", "Proportional divergence"), bty="n", pt.cex = 1.3)
legend("topleft", pch=c(0, 0), col=c("black", "black"), legend=c("", ""), bty="n", pt.cex=1.3)
axis(2, at=c(-2:3), signif(10^(-2:3), 1), las=1)
xt3 = c(1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
x3at = log10(100*log(xt3)^2)
x3atNEW = log10(100*log(xt3)^2/(2/pi))

axis(4, at=x3atNEW, signif(xt3, 4), las=1)

# General boxplots ####
ddf$logd = log10(ddf$d*100)
ddf$logd[which(ddf$logd==-Inf)]=log10(0.003907)
ddf$logd[which(ddf$logd<log10(0.003907))]=log10(0.003907) # To plot very low values at "<1.005"

tpos=-.25
xt3 = c(1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
x3at = log10(100*log(xt3)^2)
x3atNEW = log10(100*log(xt3)^2/(2/pi))

veg = ddf[ddf$tg1=="vegetative",]

medians = tapply(veg$logd, veg$tg2, median, na.rm=T)
veg$tg2 = factor(veg$tg2, levels=names(sort(medians, decreasing=T)))
levels(veg$tg2)
a = 1:9

#x11(height=5.5, width=6.5)
par(mar=c(8,4,2,2))
plot(veg$tg2, veg$logd, cex=.8, notch=F, xaxt="n", col="darkgreen", 
     ylab="Proportional divergence", xlab="", yaxt="n", ylim=c(-3,3), xlim=c(0,30))
#axis(2, at=c(-3:3), c("<0.001",c(signif(10^(c(-2:3)),1))), las=1)
axis(2, at=x3atNEW, c("<1.005", signif(xt3, 4)[-1]), las=1)

axis(1, at=a, labels = FALSE)
labels = paste0(toupper(as.character(levels(veg$tg2)))," (",tapply(veg$d>-100, veg$tg2, sum, na.rm=T),") ")
labels = sub("_"," ",labels)
text(a, par("usr")[3] + tpos, srt = 45, adj = 1,cex=.7,
     labels = labels, xpd = TRUE)
h = median(ddf$logd[ddf$tg1=="vegetative"], na.rm=T)
segments(min(a), h, max(a), h, lwd=3)

#

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

# Simpler boxplot for main text ####
#x11(height=5.5, width=5.5)
#cairo_pdf("pubfigs/FloVegBoxes.pdf", height=5.5, width=5.5, fam = "Times")
par(mar=c(8,4,2,2))
plot(veg$tg2, veg$logd, cex=.8, notch=F, xaxt="n", col="darkgreen", 
     ylab="Proportional divergence", xlab="", yaxt="n", ylim=c(-3,3), xlim=c(0,23))
axis(2, at=x3atNEW, c("<1.005", signif(xt3, 4)[-1]), las=1)

a = 1:8

axis(1, at=a, labels = FALSE)
labels = paste0(toupper(as.character(levels(veg$tg2)))," (",tapply(veg$d>-100, veg$tg2, sum, na.rm=T),") ")
labels = sub("_"," ",labels)
text(a, par("usr")[3] + tpos, srt = 45, adj = 1,cex=.7,
     labels = labels, xpd = TRUE)
h = median(ddf$logd[ddf$tg1=="vegetative"], na.rm=T)
segments(min(a), h, max(a), h, lwd=3)

# Floral
floral = ddf[ddf$tg1=="floral" & ddf$tg2!="fitness",]
#floral$logd[which(floral$logd==-Inf)]=log10(0.01)

medians = tapply(floral$logd, floral$tg2, median, na.rm=T)
floral$tg2 = factor(floral$tg2, levels=names(sort(medians,decreasing=T)))
levels(floral$tg2)
a = (max(a)+2):((max(a)+1)+length(levels(floral$tg2)))

par(new=T)
plot(floral$tg2, floral$logd, cex=.8, yaxt="n", at=a, xaxt="n", col="darkblue",
     ylab="", xlab="", ylim=c(-3,3), xlim=c(0,23))
axis(1, at=a, labels = FALSE)
labels = paste0(toupper(as.character(levels(floral$tg2)))," (",tapply(floral$logd>-100, floral$tg2, sum, na.rm=T),") ")
labels = sub("_"," ",labels)
labels[11] = "FLOWER SIZE (79)"
text(a, par("usr")[3] + tpos, srt = 45, adj = 1, cex=.7,
     labels = labels, xpd = TRUE)
h = median(ddf$logd[ddf$tg1=="floral"], na.rm=T)
segments(min(a), h, max(a), h, lwd=3)

legend("topright", pch=15, col=c("darkblue", "darkgreen"), cex=1, pt.cex=1.5, 
       bty="n", legend=c(expression(paste("Floral (", italic(n)," = 273)")), 
                         expression(paste("Vegetative (", italic(n)," = 80)"))))

#dev.off()

#### Subset data for d vs. e meta-analysis ####
ddf = na.omit(ddf)
ddf = ddf[ddf$d>0,]
ddf = ddf[ddf$evals>0,]
head(ddf)

# Check (not needed anymore)
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

#### Summary stats for univariate meta-analysis ####
length(unique(ddf$study_ID))
length(unique(ddf$species))

tapply(ddf$d, list(ddf$ms, ddf$tg1), median)

#### Univariate meta-analysis ####

# Remove some of the repeated D. scandens studies?
#ddf = ddf[ddf$study_ID!="Hansen_et_al._2003_Dalechampia_scandens_A_greenhouse",]
#ddf = ddf[ddf$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_field",]
#ddf=ddf[ddf$study_ID!="Opedal_et_al._Costa_Rica_Dalechampia_scandens_A_greenhouse",]

ddf$scale_npop = scale(log(ddf$npop), scale=F)
ddf$scale_maxdist = scale(log(ddf$maxdist), scale=F)
ddf$scale_log_evol = scale(log(ddf$evals), scale=F)

wgts = 1/((ddf$d_se^2)/(ddf$d^2)) #Mean-scaled sampling variance because analysis is on log

# Baseline meta-analytic model

# Fit with REML=T to get parameter estimates
m0 = glmmTMB(log(d)~ scale_log_evol + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m0)
exp(sqrt((2/pi)*exp(summary(m0)$coef$cond[1,1]))) #dP at intercept (mean)
r.squaredGLMM(m0)

# Refit with REML=F for AIC comparisons
m0 = glmmTMB(log(d)~ scale_log_evol + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)

# Alternative models

# Evolvability*distance interaction
m = glmmTMB(log(d)~ scale_log_evol*scale_maxdist + scale_npop + (scale_log_evol|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1]
summary(m)

# Evolvability*npop interaction
m = glmmTMB(log(d)~ scale_log_evol*scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1]
summary(m)

# Study environments
m = glmmTMB(log(d)~ scale_log_evol*environment + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1]

# Refit with REML=T to get parameter estimates
m = glmmTMB(log(d)~ scale_log_evol*environment + scale_npop + log(maxdist) + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m)
r.squaredGLMM(m)

# Suppress global intercept to get parameter estimates per group
m = glmmTMB(log(d)~ -1 + environment + scale_log_evol:environment + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m)

# Mating systems
m = glmmTMB(log(d)~ scale_log_evol*ms + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)
AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1]

# Refit with REML=T to get parameter estimates
m = glmmTMB(log(d)~ scale_log_evol*ms + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m)
r.squaredGLMM(m)

# Suppress global intercept to get parameter estimates per group
m = glmmTMB(log(d)~ -1 + dimension + scale_log_evol:dimension + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m)

# Trait dimensions (for Appendix)
m = glmmTMB(log(d)~ scale_log_evol*dimension + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights=wgts, data=ddf)

AIC(m0, m)
AIC(m0, m)$AIC[2]-AIC(m0, m)$AIC[1]

# Refit with REML=T to get parameter estimates
m = glmmTMB(log(d)~ scale_log_evol*dimension + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m)
r.squaredGLMM(m)

# Suppress global intercept to get parameter estimates per group
m = glmmTMB(log(d)~ -1 + dimension + scale_log_evol:dimension + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
#m = glmmTMB(log(d)~ scale_log_evol:dimension + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights=wgts, data=ddf)
summary(m)

# Trait groups (floral vs. vegetative)

# Subset to floral and vegetative traits only
floveg = ddf[ddf$tg1=="floral" | ddf$tg1=="vegetative",]
floveg$tg1 = factor(floveg$tg1)
floveg_wgts = 1/((floveg$d_se^2)/(floveg$d^2)) #Mean-scaled sampling variance because analysis is on log

# New baseline model for AIC comparison
m0 = glmmTMB(log(d) ~ scale_log_evol + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights = floveg_wgts, data=floveg)

m = glmmTMB(log(d) ~ scale_log_evol*tg1 + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights = floveg_wgts, data=floveg)
m1 = glmmTMB(log(d) ~ scale_log_evol+tg1 + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=F, weights = floveg_wgts, data=floveg)

AIC(m0, m, m1)
AIC(m0, m)$AIC[2] - AIC(m0, m)$AIC[1]

# Refit with REML=T to get parameter estimates
m = glmmTMB(log(d) ~ scale_log_evol*tg1 + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights = floveg_wgts, data=floveg)
summary(m)
r.squaredGLMM(m)

# Suppress global intercept to get parameter estimates per group
m = glmmTMB(log(d) ~ -1 + tg1 + scale_log_evol:tg1 + scale_npop + scale_maxdist + (scale_log_evol|study_ID) + (1|species), REML=T, weights = floveg_wgts, data=floveg)
summary(m)

#### Plotting evolvability vs. divergence ####

# All data
par(mfrow=c(1,1))
plot(log10(ddf$evals),ddf$logd,
     xlab="Evolvability (%)",
     ylab="Among-population variance (%)",
     pch=1,cex=1*sqrt(ddf$npop),
     col="black",
     xlim=c(-2.5, 2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))
axis(2,c(-4,-3,-2,-1,0,1,2),10^c(-4,-3,-2,-1,0,1,2), las=1)

quantile(ddf$d)
points(tapply(ddf$logd, ddf$study_ID, mean, na.rm=T),tapply(log10(ddf$evals), ddf$study_ID, mean, na.rm=T), pch=16)
summary(lm(tapply(ddf$logd, ddf$study_ID, mean, na.rm=T)~tapply(log(ddf$evals), ddf$study_ID, mean, na.rm=T)))

# Floral and vegetative
floveg = ddf[ddf$tg1=="floral" | ddf$tg1=="vegetative",]
floveg$tg1 = factor(floveg$tg1)
m = glmmTMB(log(d*100) ~ -1 + tg1 + log(evals):tg1 + scale_npop+ scale_maxdist 
         + (log(evals)|study_ID) +(1|species), weights = floveg_wgts, data=floveg)

xt3 = c(1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
#x3at = log10(100*log(xt3)^2)
x3at = log10(100*log(xt3)^2/(2/pi))

cols = c(rgb(0, 0, 0.545, .35), 
         rgb(0.004, 0.196, 0.125, .35),
         rgb(.545, 0, 0, .35))

cairo_pdf("pubfigs/univar_evol_div.pdf", width=4, height=8, fam="Times")
#x11(width=4, height=8)
par(mfrow=c(3,1), mar=c(4,5,0.1,2))

plot(log10(floveg$evals), floveg$logd,
     xlab="",
     ylab="",
     pch=1, cex=1*sqrt(floveg$npop),
     col=cols[1:2][as.numeric(floveg$tg1)],
     xlim=c(-2.5, 3.5), ylim=c(-3,3), xaxt="n", yaxt="n")
axis(1, c(-2,-1,0,1,2), 10^c(-2,-1,0,1,2))
#axis(2, c(-3,-2,-1,0,1,2), c("<0.001", 10^c(-2,-1,0,1,2)), las=1)
axis(2, at=x3at, c("<1.005", signif(xt3, 4)[-1]), las=1)
mtext("Proportional divergence", 2, line=3.3, cex=0.8)

par(new=T)
plot(floveg$tg1, floveg$logd, at=c(2.5, 3), boxwex=0.4, xlab="", ylab="",
     col=cols[1:2],
     xlim=c(-2.5, 3.5), ylim=c(-3,3), xaxt="n", yaxt="n")

x1=seq(min(log(floveg$evals[floveg$tg1=="floral"])),
       max(log(floveg$evals[floveg$tg1=="floral"])), .1)
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[5,1]*x1 + 
  summary(m)$coef$cond[3,1]*mean(ddf$scale_npop[ddf$tg1=="floral"]) +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="floral"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e = log(median(ddf$evals[ddf$tg1=="floral"]))
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[5,1]*me_e +
  summary(m)$coef$cond[3,1]*mean(ddf$scale_npop[ddf$tg1=="floral"]) +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="floral"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(floveg$evals[floveg$tg1=="vegetative"])),
       max(log(floveg$evals[floveg$tg1=="vegetative"])), .1)
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[6,1]*x2 +
  summary(m)$coef$cond[3,1]*mean(ddf$scale_npop[ddf$tg1=="vegetative"]) +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="vegetative"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$tg1=="vegetative"]))
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[6,1]*me_e +
  summary(m)$coef$cond[3,1]*mean(ddf$scale_npop[ddf$tg1=="vegetative"]) +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_maxdist[ddf$tg1=="vegetative"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

legend("topleft", pch=15, col=c("darkblue", "darkgreen"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Floral", "Vegetative"))

# Mating systems
m = glmmTMB(log(d*100) ~ -1 + ms + log(evals):ms + scale_npop+ scale_maxdist 
         + (log(evals)|study_ID) +(1|species), weight = wgts, data=ddf)
summary(m)$coef$cond

plot(log10(ddf$evals), ddf$logd,
     xlab="",
     ylab="",
     pch=1, cex=1*sqrt(ddf$npop),
     col=cols[as.numeric(ddf$ms)],
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")
axis(1, c(-2,-1,0,1,2), 10^c(-2,-1,0,1,2))
#axis(2, c(-3,-2,-1,0,1,2), c("<0.001", 10^c(-2,-1,0,1,2)), las=1)
axis(2, at=x3at, c("<1.005", signif(xt3, 4)[-1]), las=1)
#mtext("Evolvability (%)", 1, line=2.5)
mtext("Proportional divergence", 2, line=3.3, cex=0.8)

par(new=T)
plot(ddf$ms, ddf$logd, at=c(2.5, 3, 3.5), boxwex=0.4, xlab="", ylab="",
     col=c(rgb(0, 0, 0.545, .35), rgb(0.004, 0.196, 0.125, .35), rgb(.545, 0, 0, .35)),
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")

x1=seq(min(log(ddf$evals[ddf$ms=="S"])),
       max(log(ddf$evals[ddf$ms=="S"])), .1)
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[6,1]*x1 +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$ms=="S"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$ms=="S"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e=log(median(ddf$evals[ddf$ms=="S"]))
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[6,1]*me_e +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$ms=="S"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$ms=="S"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(ddf$evals[ddf$ms=="M"])),
       max(log(ddf$evals[ddf$ms=="M"])), .1)
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[7,1]*x2 +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$ms=="M"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$ms=="M"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$ms=="M"])) 
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[7,1]*me_e +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$ms=="M"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$ms=="M"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

x3=seq(min(log(ddf$evals[ddf$ms=="O"])),
       max(log(ddf$evals[ddf$ms=="O"])), .1)
y = summary(m)$coef$cond[3,1] + summary(m)$coef$cond[8,1]*x3 +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$ms=="O"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$ms=="O"])
lines(log10(exp(x3)), log10(exp(y)), col ="darkred", lwd=2)

me_e = log(median(ddf$evals[ddf$ms=="O"]))
y = summary(m)$coef$cond[3,1] + summary(m)$coef$cond[8,1]*me_e +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$ms=="O"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$ms=="O"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkred")

legend("topleft", pch=15, col=c("darkblue", "darkgreen", "darkred"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Selfing", "Mixed", "Outcrossing"))

# Environments
m = glmmTMB(log(d*100) ~ -1 + environment + log(evals):environment + scale_npop+ scale_maxdist 
         + (log(evals)|study_ID) +(1|species), weights=wgts, data=ddf)

plot(log10(ddf$evals), ddf$logd,
     xlab="",
     ylab="",
     pch=1, cex=1*sqrt(ddf$npop),
     col=cols[as.numeric(ddf$environment)],
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")
axis(1, c(-2,-1,0,1,2), 10^c(-2,-1,0,1,2))
axis(2, at=x3at, c("<1.005", signif(xt3, 4)[-1]), las=1)
mtext("Evolvability (%)", 1, line=2.5, cex=0.8)
mtext("Proportional divergence", 2, line=3.3, cex=0.8)

par(new=T)
plot(ddf$environment, ddf$logd, at=c(2.5, 3, 3.5), boxwex=0.4, xlab="", ylab="",
     col=c(rgb(0, 0, 0.545, .35), rgb(0.004, 0.196, 0.125, .35), rgb(.545, 0, 0, .35)),
     xlim=c(-2.5, 3.75), ylim=c(-3.5,3), xaxt="n", yaxt="n")

x1=seq(min(log(ddf$evals[ddf$environment=="greenhouse"])),
       max(log(ddf$evals[ddf$environment=="greenhouse"])), .1)
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[6,1]*x1 +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$environment=="greenhouse"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$environment=="greenhouse"])
lines(log10(exp(x1)), log10(exp(y)), col ="darkblue", lwd=2)

me_e=log(median(ddf$evals[ddf$environment=="greenhouse"]))
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[6,1]*me_e +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$environment=="greenhouse"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$environment=="greenhouse"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkblue")

x2=seq(min(log(ddf$evals[ddf$environment=="common_garden"])),
       max(log(ddf$evals[ddf$environment=="common_garden"])), .1)
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[7,1]*x2 +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$environment=="common_garden"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$environment=="common_garden"])
lines(log10(exp(x2)), log10(exp(y)), col ="darkgreen", lwd=2)

me_e = log(median(ddf$evals[ddf$environment=="common_garden"])) 
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[7,1]*me_e +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$environment=="common_garden"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$environment=="common_garden"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkgreen")

x3=seq(min(log(ddf$evals[ddf$environment=="field"])),
       max(log(ddf$evals[ddf$environment=="field"])), .1)
y = summary(m)$coef$cond[3,1] + summary(m)$coef$cond[8,1]*x3 +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$environment=="field"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$environment=="field"])
lines(log10(exp(x3)), log10(exp(y)), col ="darkred", lwd=2)

me_e = log(median(ddf$evals[ddf$environment=="field"]))
y = summary(m)$coef$cond[3,1] + summary(m)$coef$cond[8,1]*me_e +
  summary(m)$coef$cond[4,1]*mean(ddf$scale_npop[ddf$environment=="field"]) +
  summary(m)$coef$cond[5,1]*mean(ddf$scale_maxdist[ddf$environment=="field"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="darkred")

legend("topleft", pch=15, col=c("darkblue", "darkgreen", "darkred"), cex=1, pt.cex=1.5, 
       bty="n", legend=c("Greenhouse", "Garden", "Field"))

dev.off()

# Dimensions
tapply(ddf$dimension, ddf$dimension, length)

m = glmmTMB(log(d*100) ~ -1 + dimension + log(evals):dimension + scale_npop + scale_maxdist 
         + (log(evals)|study_ID) +(1|species), weights = wgts, data=ddf)

cols = c("blue3", "firebrick", "green3", "yellow3","cyan", "black") 

x11(height=5, width=6)
par(mfrow=c(1,1), mar=c(4,4,2,7))
plot(log10(ddf$evals), ddf$logd,
     xlab="Evolvability (%)",
     ylab="Proportional divergence",
     pch=1, cex=1*sqrt(ddf$npop),
     col=cols[as.numeric(ddf$dimension)],
     xlim=c(-2.5, 2),ylim=c(-4,3),xaxt="n", yaxt="n")
axis(1,c(-2,-1,0,1,2,3),10^c(-2,-1,0,1,2,3))

xt3 = c(1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
x3at = log10(100*log(xt3)^2/(2/pi))
axis(2, at=x3at, c("<1.005", signif(xt3, 4)[-1]), las=1)


x1=seq(min(log(ddf$evals[ddf$dimension=="ratio"])),
       max(log(ddf$evals[ddf$dimension=="ratio"])), .1)
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[9,1]*x1 + 
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="ratio"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="ratio"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[1], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="linear"])),
       max(log(ddf$evals[ddf$dimension=="linear"])), .1)
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[10,1]*x1 + 
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="linear"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="linear"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[2], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="time"])),
       max(log(ddf$evals[ddf$dimension=="time"])), .1)
y = summary(m)$coef$cond[3,1] + summary(m)$coef$cond[11,1]*x1 + 
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="time"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="time"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[3], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="mass_volume"])),
       max(log(ddf$evals[ddf$dimension=="mass_volume"])), .1)
y = summary(m)$coef$cond[4,1] + summary(m)$coef$cond[12,1]*x1 + 
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="mass_volume"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="mass_volume"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[4], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="area"])),
       max(log(ddf$evals[ddf$dimension=="area"])), .1)
y = summary(m)$coef$cond[5,1] + summary(m)$coef$cond[13,1]*x1 + 
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="area"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="area"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[5], lwd=2)

x1=seq(min(log(ddf$evals[ddf$dimension=="count"])),
       max(log(ddf$evals[ddf$dimension=="count"])), .1)
y = summary(m)$coef$cond[6,1] + summary(m)$coef$cond[14,1]*x1 + 
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="count"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="count"])
lines(log10(exp(x1)), log10(exp(y)), col=cols[6], lwd=2)

me_e = log(median(ddf$evals[ddf$dimension=="ratio"]))
y = summary(m)$coef$cond[1,1] + summary(m)$coef$cond[9,1]*me_e +
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="ratio"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="ratio"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="linear"]))
y = summary(m)$coef$cond[2,1] + summary(m)$coef$cond[10,1]*me_e +
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="linear"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="linear"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="time"]))
y = summary(m)$coef$cond[3,1] + summary(m)$coef$cond[11,1]*me_e +
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="time"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="time"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="mass_volume"]))
y = summary(m)$coef$cond[4,1] + summary(m)$coef$cond[12,1]*me_e +
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="mass_volume"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="mass_volume"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="area"]))
y = summary(m)$coef$cond[5,1] + summary(m)$coef$cond[13,1]*me_e +
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="area"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="area"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

me_e = log(median(ddf$evals[ddf$dimension=="count"]))
y = summary(m)$coef$cond[6,1] + summary(m)$coef$cond[14,1]*me_e +
  summary(m)$coef$cond[7,1]*mean(ddf$scale_npop[ddf$dimension=="count"]) +
  summary(m)$coef$cond[8,1]*mean(ddf$scale_maxdist[ddf$dimension=="count"])
points(log10(exp(me_e)), log10(exp(y)), pch=16, col="black", cex=1.3)

legend(x=2.2, y = 2, pch=15, col=cols, cex=1, pt.cex=1.5, 
       bty="n", legend=levels(ddf$dimension), xpd=T)

# Removing the dimension effect by taking residuals
resd = summary(lm(log10(ddf$d*100)~ddf$dimension))$residuals
rese = summary(lm(log10(ddf$evals)~ddf$dimension))$residuals

plot(rese, resd)
m = glmmTMB(resd ~ rese + scale_npop + scale_maxdist
         + (rese|study_ID) + (1|species), weights = wgts, data=ddf)
summary(m)

r.squaredGLMM(m)

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


#### END OF ANALYSES REPORTED IN PAPER ####

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
plotSubset("dimension", "time")

# Individual studies
studies = unique(ddf$study_ID)
pdf("figs/univariate_G_D_plots2.pdf")
for(s in 1:length(studies)){
  study=studies[s]
  plotSubset("study_ID", study)
}
dev.off()

#### Test off-setting of area and cubic measures
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


# Outcrossing rates ####
for(i in 1:nrow(ddf)){
  if(!is.na(ddf$ms[i])){
    if(is.na(ddf$tmvals[i]) & ddf$ms[i]=="O")
      ddf$tmvals[i]=1
  }}

plot(ddf$tmvals, log(ddf$d))

#Scaling of trait variances with the mean
#ddat=ddat[ddat$trait!="organ_size",]
#plot(log10(ddat$mean), log10(ddat$sd))
#lines(-10:10, -10:10)


flo = ddf[which(ddf$tg1=="floral"),]
veg = ddf[which(ddf$tg1=="vegetative"),]
fit = ddf[which(ddf$tg2=="fit"),]
size = ddf[which(ddf$tg2=="flowersize"),]

length(which(fit$evals<median(ddf$evals, na.rm=T)))/sum(fit$evals>-Inf, na.rm=T)
length(which(fit$d<median(ddf$d, na.rm=T)))/sum(fit$d>-Inf, na.rm=T)
length(which(size$evals<median(ddf$evals, na.rm=T)))/sum(size$evals>-Inf, na.rm=T)
length(which(size$d<median(ddf$d, na.rm=T)))/sum(size$d>-Inf, na.rm=T)

par(mfrow=c(1,1))
plot(log10(flo$evals), log10(flo$d), ylim=c(-7.2, 0.5), yaxt="n")

plot(log10(ddf$evals), log10(ddf$d), ylim=c(-6, 1), col="lightgrey",
     xlab="Evolvability", ylab="Proportional divergence", yaxt="n")
points(log10(fit$evals), log10(fit$d), pch=16)
abline(h=log10(median(fit$d)))
abline(v=log10(median(fit$evals, na.rm=T)))

xt3 = c(1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 3)
x3at = log10(log(xt3)^2)
axis(2, at=x3at, signif(xt3, 4), las=1)


# Back-transformation
dp = exp(sqrt(ddf$d*(2/pi)))
cor(ddf$d, (log(dp))^2/(2/pi))

# Contextual model
emeans = tapply(log(ddf$evals), ddf$study_ID, mean)

meanvals=NULL
for(i in 1:nrow(ddf)){
  meanvals[i]=emeans[which(names(emeans)==ddf$study_ID[i])] 
}
data.frame(meanvals, ddf$study_ID)

m = lmer(log(d*100) ~ log(evals) + meanvals + scale_maxdist + scale_npop +
           + (1|species/study_ID), data=ddf)
summary(m)

# Formal meta-analysis using Almer_SE
SE = sqrt((ddf$d_se^2)/(ddf$d^2))
plot(ddf$npop, SE)
m0 = evolvability::Almer_SE(log(d) ~ log(evals) + log(npop) + log(maxdist) + (log(evals)|study_ID) + (1|species), 
                            SE=SE, maxiter = 100, data=ddf)
summary(m0)

m1 = evolvability::Almer_SE(log(d) ~ log(evals) + log(npop) + log(maxdist) + (log(evals)|study_ID), 
                            SE=SE, maxiter = 100, data=ddf)
AIC(m0, m1)

#### Formal meta-analysis using MCMCglmm

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


