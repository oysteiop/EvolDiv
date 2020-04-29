################################################################################
##### Combine the results of the eigenvector and divergence vector analyses ####
################################################################################
rm(list=ls())
library(plyr)
library(lme4)

load(file="gdDat.RData")
load(file="deltaDat.RData")
load("data/EVOBASE.RData")
load("data/POPBASE.RData")

names(gdDat)
gdMeans= ddply(gdDat, .(species), summarize,
                        emean = median(emean),
                        emax = median(emax),
                        cmean = median(cmean),
                        dmean= median(dmean),
                        betaG = median(betaG),
                        r2G = median(r2G),
                        betaD = median(betaD),
                        r2D = median(r2D),
                        thetaGD = median(theta),
                        r2All = median(r2All),
                        ntraits = median(ntraits))

deltaMeans = ddply(deltaDat, .(species), summarize,
                emean = median(emean),
                cmean = median(cmean),
                edelta= median(edelta),
                cdelta = median(cdelta),
                div = median(div),
                theta_delta=median(theta))

comb = merge(gdMeans, deltaMeans, by="species", all=T)

gdDat$ID = paste(gdDat$g, gdDat$d, sep="_")
deltaDat$ID = paste(deltaDat$g, deltaDat$d, sep="_")

gdMeans2= ddply(gdDat, .(species, ID), summarize,
               emean = median(emean),
               emax = median(emax),
               cmean = median(cmean),
               dmean= mean(dmean),
               betaG = median(betaG),
               betaG_SE = betaG_SE,
               r2G = median(r2G),
               betaD = median(betaD),
               r2D = median(r2D),
               thetaGD = median(theta),
               r2All = median(r2All),
               ntraits = median(ntraits),
               dims = dims,
               ndims = ndims,
               traitgroups = traitgroups)

deltaMeans2 = ddply(deltaDat, .(species, ID), summarize,
                   g = unique(g)[1],
                   d = unique(d)[1],
                   emean = median(emean),
                   cmean = median(cmean),
                   edelta= median(edelta),
                   cdelta = median(cdelta),
                   div = mean(div),
                   theta_delta = median(theta),
                   dims = unique(dims),
                   ndims = mean(ndims),
                   traitgroups = unique(traitgroups))

comb2 = merge(gdMeans2, deltaMeans2, by="ID", all=T)
head(comb2, 5)

comb2$ndims = apply(subset(comb2, select=c("ndims.x", "ndims.y")), 1, mean, na.rm=T) 
comb2$dims = unlist(apply(subset(comb2, select=c("dims.x", "dims.y")), 1, function(x) sort(unique(x))))
comb2$traitgroups = unlist(apply(subset(comb2, select=c("traitgroups.x", "traitgroups.y")), 1, function(x) sort(unique(x))))

table(comb2$traitgroups)
table(comb2$dims)
table(comb2$ndims)

x11(height=4.5, width=11)
par(mfrow=c(1,3))
par(mar=c(8,4,2,2))

plot(factor(comb2$ndims), comb2$betaG, xlab="Number of trait dimensions", ylab="Slope for G eigenvectors", las=1)
abline(h=1, lty=2)

medians = tapply(comb2$betaG, comb2$dims, median, na.rm=T)
comb2$dims = factor(comb2$dims, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$dims), comb2$betaG, xaxt="n", xlab="", las=1,
     ylab="Slope for G eigenvectors")
axis(1, at=1:10, labels=rep("", 10))
text(1:10, par("usr")[3] - .05, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$dims), xpd = TRUE)

medians = tapply(comb2$betaG, comb2$traitgroups, median, na.rm=T)
comb2$traitgroups = factor(comb2$traitgroups, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$traitgroups), comb2$betaG, xaxt="n", xlab="", las=1,
     ylab="Slope for G eigenvectors")
axis(1, at=1:4, labels=rep("", 4))
text(1:4, par("usr")[3] - .05, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$traitgroups), xpd = TRUE)


x11(height=4.5, width=11)
par(mfrow=c(1,3))
par(mar=c(8,4,2,2))

plot(factor(comb2$ndims), comb2$edelta/comb2$emean.y, xlab="Number of trait dimensions", 
     ylab="edelta/emean", las=1)
abline(h=1, lty=2)

medians = tapply(comb2$edelta/comb2$emean.y, comb2$dims, median, na.rm=T)
comb2$dims = factor(comb2$dims, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$dims), comb2$edelta/comb2$emean.y, xaxt="n", xlab="", las=1,
     ylab="edelta/emean")
axis(1, at=1:10, labels=rep("", 10))
text(1:10, par("usr")[3] - .05, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$dims), xpd = TRUE)

medians = tapply(comb2$edelta/comb2$emean.y, comb2$traitgroups, median, na.rm=T)
comb2$traitgroups = factor(comb2$traitgroups, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$traitgroups), comb2$edelta/comb2$emean.y, xaxt="n", xlab="", las=1,
     ylab="edelta/emean")
axis(1, at=1:4, labels=rep("", 4))
text(1:4, par("usr")[3] - .05, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$traitgroups), xpd = TRUE)

#Mating systems from EVOBASE
ms = NULL
for(i in 1:nrow(comb2)){
ms[i] = paste0(EVOBASE[[as.character(comb2$g[i])]]$MS, collapse = "+")
}
ms = factor(ms, levels=c("S", "M", "O"))
table(ms)
data.frame(comb2$ID, ms)

comb2$ms = ms
table(comb2$ms)

plot(ms, log10(comb2$div))
plot(ms, comb2$betaG)
plot(ms, comb2$edelta/comb2$emean.y)


m = lmer(betaG ~ ms + (1|species.y), weights=(1/comb2$betaG_SE^2),data=comb2)
summary(m)
m2 = Almer_SE(betaG ~ ms + (1|species.y), SE=comb2$betaG_SE, data=comb2)
summary(m2)

# Study environments (for the divergence data) from POPBASE
env = NULL
for(i in 1:nrow(comb2)){
  env[[i]] = paste0(POPBASE[[as.character(comb2$d[i])]]$Env, collapse="+")
}
env = factor(env, levels=c("field", "common_garden", "greenhouse"))

table(env)
data.frame(comb2$ID, env)

comb2$env = env
table(comb2$env)

plot(comb2$env, log10(comb2$dmean))
plot(comb2$env, comb2$betaG)
plot(comb2$env, log(comb2$edelta/comb2$emean.y))
abline(h=0)

m = lmer(betaG ~ env + (1|species.y), weights=(1/comb2$betaG_SE^2), REML=T, data=comb2)
summary(m)

m2 = Almer_SE(betaG ~ env + (1|species.y), SE=comb2$betaG_SE, data=comb2)
summary(m2)

#### Plots comparing GD and delta analyses ####
x11(height=4, width=7)
par(mfrow=c(1,2))
plot(log(comb2$edelta/comb2$emean.x), comb2$betaG, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="Slope of G eigenvectors")
points(log(comb$edelta/comb$emean.x), comb$betaG, pch=16)
abline(v=0, lty=2)
abline(h=1, lty=2)

plot(log(comb2$cdelta/comb2$cmean.x), comb2$betaG, pch=16, col="lightgrey", las=1,
     xlab="log(c_delta/c_mean)", ylab="Slope of G eigenvectors")
points(log(comb$cdelta/comb$cmean.x), comb$betaG, pch=16)
abline(v=0, lty=2)
abline(h=1, lty=2)


plot(log(comb2$edelta/comb2$emean.x), comb2$r2G, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="R2 for G eigenvectors")
points(log(comb$edelta/comb$emean.x), comb$r2G, pch=16)
abline(v=0, lty=2)


par(mfrow=c(1,1))
plot(log(comb2$edelta/comb2$emean.x), comb2$thetaGD, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="Angle between gmax and dmax")
points(log(comb$edelta/comb$emean.x), comb$thetaGD, pch=16)
abline(v=0, lty=2)

plot(comb2$edelta/comb2$emax, comb2$thetaGD, pch=16, col="lightgrey", las=1,
     xlab="e_delta/e_max", ylab="Angle between gmax and dmax")
points(comb$edelta/comb$emax, comb$thetaGD, pch=16)
#abline(v=0, lty=2)

plot(comb2$r2All, comb2$thetaGD, pch=16, col="lightgrey", las=1, xlim=c(0,1), ylim=c(0,90),
     xlab="Overall R2", ylab="Angle between gmax and dmax")
points(comb$r2All, comb$thetaGD, pch=16)
#abline(v=0, lty=2)

plot(comb2$r2All, comb2$thetaGD, pch=16)


plot(log(comb2$cdelta/comb2$cmean.x), comb2$thetaGD, pch=16, col="lightgrey", las=1,
     xlab="log(c_delta/c_mean)", ylab="Angle between gmax and dmax")
points(log(comb$cdelta/comb$cmean.x), comb$thetaGD, pch=16)
abline(v=0, lty=2)



plot(log(comb2$edelta/comb2$emean.x), comb2$r2All, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="Angle between gmax and dmax")
points(log(comb$edelta/comb$emean.x), comb$r2All, pch=16)
abline(v=0, lty=2)


cor.test(log(comb$edelta/comb$emean.x), comb$betaG)
cor.test(log(comb2$edelta/comb2$emean.x), comb2$betaG)

cor.test(log(comb$cdelta/comb$cmean.x), comb$betaG)
cor.test(log(comb2$cdelta/comb2$cmean.x), comb2$betaG)

cor.test(log(comb$edelta/comb$emean.x), comb$thetaGD)
cor.test(log(comb2$edelta/comb2$emean.x), comb2$thetaGD)

summary(lmer(betaG ~ log(edelta/emean.x) + (1|sp.x), data=comb2))
summary(lmer(betaG ~ log(cdelta/cmean.x) + (1|sp.x), data=comb2))
summary(lmer(thetaGD ~ log(edelta/emean.x) + (1|sp.x), data=comb2))
summary(lmer(thetaGD ~ log(cdelta/cmean.x) + (1|sp.x), data=comb2))


x11(height=7, width=7)
par(mfrow=c(2,2))
plot(log(comb2$edelta/comb2$emean.x), comb2$betaD, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="Slope of D eigenvectors")
points(log(comb$edelta/comb$emean.x), comb$betaD, pch=16)

plot(log(comb2$edelta/comb2$emean.x), comb2$r2D, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="R2 for D eigenvectors")
points(log(comb$edelta/comb$emean.x), comb$r2D, pch=16)

plot(log(comb2$cdelta/comb2$cmean.x), comb2$betaD, pch=16, col="lightgrey", las=1,
     xlab="log(c_delta/c_mean)", ylab="Slope of D eigenvectors")
points(log(comb$cdelta/comb$cmean.x), comb$betaD, pch=16)

plot(log(comb2$edelta/comb2$emean.x), comb2$thetaGD, pch=16, col="lightgrey", las=1,
     xlab="log(e_delta/e_mean)", ylab="Angle between gmax and dmax")
points(log(comb$edelta/comb$emean.x), comb$thetaGD, pch=16)


