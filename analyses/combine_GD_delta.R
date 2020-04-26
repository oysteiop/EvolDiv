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
               r2G = median(r2G),
               betaD = median(betaD),
               r2D = median(r2D),
               thetaGD = median(theta),
               r2All = median(r2All),
               ntraits = median(ntraits),
               traitgroup = unique(traitgroup)[1])

deltaMeans2 = ddply(deltaDat, .(species, ID), summarize,
                   g = unique(g)[1],
                   d = unique(d)[1],
                   emean = median(emean),
                   cmean = median(cmean),
                   edelta= median(edelta),
                   cdelta = median(cdelta),
                   div = mean(div),
                   theta_delta = median(theta))

comb2 = merge(gdMeans2, deltaMeans2, by="ID", all=T)
head(comb2)

ms = NULL
for(i in 1:nrow(comb2)){
ms[i] = EVOBASE[[comb2$g[i]]]$MS
}
ms = factor(ms, levels=c("S", "M", "O"))
plot(ms, comb2$betaG)
plot(ms, comb2$betaD)
plot(ms, comb2$r2all)

m = lmer((edelta/emean.y)~ms +(1|species.y), data=comb2)
summary(m)

tg = list()
for(i in 1:nrow(comb2)){
  tg[[i]] = paste0(substr(sort(unique(EVOBASE[[comb2$g[i]]]$Groups)), 1, 3), collapse="+")
}
head(tg)

tg=unlist(tg)
table(tg)

plot(factor(tg), comb2$edelta/comb2$emean.x)

env = list()
for(i in 1:nrow(comb2)){
  env[[i]] = POPBASE[[comb2$d[i]]]$Env
}
head(env)

env=unlist(env)
table(env)

plot(factor(env), comb2$edelta/comb2$emean.x)
plot(factor(env), comb2$betaG)



comb2$traitgroup=factor(comb2$traitgroup)
par(mfrow=c(1,1))
plot(comb2$traitgroup, comb2$betaG)
plot(comb2$traitgroup, comb2$betaD)
plot(comb2$traitgroup, comb2$r2G)
plot(comb2$traitgroup, comb2$edelta/comb2$emean.x)
plot(comb2$traitgroup, comb2$betaG)

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


