################################################################################
##### Postprocessing and figures for GD-analyses ####
################################################################################

rm(list=ls())

library(lme4)

load(file="gdDatSE.RData")

# Summary stats ####
sum(gdDat$betaG>-Inf) # Number of multivariate scaling relationships for G/D directions

median(gdDat$betaD) # Median slope for D directions
median(gdDat$r2D) # Median r2 for D directions

median(gdDat$betaG) # Median slope for G directions
median(gdDat$r2G) # Median r2 for D directions

sum(gdDat$betaP>-Inf, na.rm=T) # Number of multivariate scaling relationships for P directions
median(gdDat$betaP, na.rm=T)
median(gdDat$r2P, na.rm=T)

median(gdDat$betaD_cond) # Median slope for D directions (conditional evolvability)
median(gdDat$betaP_cond, na.rm=T) # # Median slope for P directions (conditional evolvability)

# Weighed mean
wval = NULL
u = NULL
for(i in 1:nrow(gdDat)){
  u[i] = 1/(gdDat$betaG_SE[i]^2)
  wval[i] = gdDat$betaG[i]*u[i]
}
sum(wval, na.rm=T)/sum(u, na.rm=T)
mean(gdDat$betaG)

tapply(gdDat$g, gdDat$d, length)

# Meta-analysis
m = lmer(betaD~1 + (1|species/d), weights= (1/(betaD_SE^2)), data=gdDat[which(gdDat$betaD_SE>0),])
summary(m)

m = lmer(betaG~1 + (1|species/d), weights=(1/(betaG_SE^2)), data=gdDat[which(gdDat$betaG_SE>0),])
summary(m)

m = lmer(betaP~1 + (1|species/d), weights= (1/(betaP_SE^2)), data=gdDat[which(gdDat$betaP_SE>0),])
summary(m)

m = lmer(betaD_cond~1 + (1|species/d), weights=(1/(betaD_cond_SE^2)), data=gdDat[which(gdDat$betaD_cond_SE>0),])
summary(m)

m = lmer(betaP_cond~1 + (1|species/d), weights=(1/(betaP_cond_SE^2)), data=gdDat[which(gdDat$betaP_cond_SE>0),])
summary(m)

#aDat = gdDat[which(gdDat$betaG_SE>0),]
#mSE = Almer_SE(betaG ~ 1 + (1|species/d), SE=aDat$betaG_SE, data=aDat)
#m=lm(betaG~1, weights=(1/gdDat$betaG_SE^2), data=gdDat)
#summary(mSE)

# Means per species
meanDat = ddply(gdDat, .(species), summarize,
                emean = median(emean),
                cmean = median(cmean),
                dmean= median(dmean),
                betaG = median(betaG),
                r2G = median(r2G),
                betaD = median(betaD),
                r2D = median(r2D),
                betaD_cond = median(betaD_cond),
                theta = median(theta))
meanDat

min(c(range(gdDat$betaG), range(gdDat$betaD)))
max(c(range(gdDat$betaG), range(gdDat$betaD)))

#cm = cov(cbind(gdDat$betaG, gdDat$betaD))/38
#library(ellipse)
#ellipse(cm, centre=c(mean(gdDat$betaG), mean(gdDat$betaD)))

#gdDat = gdDat[gdDat$ndims==1,]

# Slope scatterplot figure ####
cairo_pdf("pubfigs/slopescatter.pdf", height=6, width=9, family="Times")
#x11(height=6, width=9)
mat = matrix(c(1,2,3,4,5,5,6,6,5,5,6,6), nrow=3, byrow=T)
layout(mat = mat)
par(mar=c(2,4,4,2), oma=c(1,0,0,0))
hist(gdDat$r2G, breaks=10, xlab="", ylab="", main=expression(paste(r^2, " G-directions")), las=1)
mtext("Frequency", 2, line=2.5)
text(-.275, 10, "(A)", cex=1.5, xpd=T)
hist(gdDat$r2D, breaks=10, xlab="", ylab="", main=expression(paste(r^2, " D-directions")), las=1)
hist(gdDat$r2P, breaks=10, xlab="", ylab="", main=expression(paste(r^2, " P-directions")), las=1)
hist(gdDat$r2T, breaks=10, xlab="", ylab="", main=expression(paste(r^2, " original traits")), las=1)

par(mar=c(4,4,1,2))
plot(gdDat$betaG, gdDat$betaD, cex=gdDat$r2All*6, lwd=2, col="lightgrey",
     ylim=c(-1,4), xlim=c(-.2,2), xlab="", ylab="", las=1)
#segments(gdDat$betaG-gdDat$betaG_SE, gdDat$betaD, gdDat$betaG+gdDat$betaG_SE, gdDat$betaD, col="grey")
#segments(gdDat$betaG, gdDat$betaD-gdDat$betaD_SE, gdDat$betaG, gdDat$betaD+gdDat$betaD_SE, col="grey")
points(meanDat$betaG, meanDat$betaD, pch=16)
points(median(meanDat$betaG), median(meanDat$betaD), pch=16, col="blue", cex=1.5)
abline(h=1, lty=2)
abline(v=1, lty=2)
mtext("Slope for G-directions", 1, line=3)
mtext("Slope for D-directions", 2, line=2.5)
#points(gdDat$betaG, gdDat$betaP, col="lightblue", pch=16)
lines(-10:10, -10:10, lty=2)
text(-.5, 4, "(B)", cex=1.5, xpd=T)

#lines(ellipse(cm, centre=c(mean(gdDat$betaG), mean(gdDat$betaD))))

plot(gdDat$r2G, gdDat$betaG, pch=16, las=1, xlim=c(0,1), ylim=c(0, 3), ylab="", xlab="")
points(gdDat$r2D, gdDat$betaD, pch=1, col="black")
points(gdDat$r2P, gdDat$betaP, pch=16, col="firebrick")
#points(gdDat$r2T, gdDat$betaT, pch=16, col="blue")

mtext(expression(paste(r^2, "")), 1, line=3)
mtext("Slope", 2, line=2.5)
legend("topleft", pch=c(16,1,16), cex=1.3, col=c("black", "black", "firebrick"), legend=c("G-directions", "D-directions", "P-directions"))
abline(h=1)
text(-.165, 3, "(C)", cex=1.5, xpd=T)

dev.off()

# Models to explain variation in betaG
library(glmmTMB)
library(MuMIn)
m1 = glmmTMB(betaG ~ r2All + ntraits + log(dmean) + (1|species), weights=1/betaG_SE^2, data=gdDat)
m2 = glmmTMB(betaG ~ r2All + ntraits + (1|species), weights=1/betaG_SE^2, data=gdDat)
m3 = glmmTMB(betaG ~ ntraits + log(dmean) + (1|species), weights=1/betaG_SE^2, data=gdDat)
m4 = glmmTMB(betaG ~ r2All + log(dmean) + (1|species), weights=1/betaG_SE^2, data=gdDat)
m5 = glmmTMB(betaG ~ r2All + (1|species), weights=1/betaG_SE^2, data=gdDat)
AIC(m1, m2, m3, m4, m5)
summary(m1)

m1 = glmmTMB(betaG ~ r2All + ntraits + log(dmean) + (1|species), data=gdDat)
m2 = glmmTMB(betaG ~ r2All + ntraits + (1|species), data=gdDat)
m3 = glmmTMB(betaG ~ ntraits + log(dmean) + (1|species), data=gdDat)
m4 = glmmTMB(betaG ~ r2All + log(dmean) + (1|species), data=gdDat)
m5 = glmmTMB(betaG ~ r2All + (1|species), data=gdDat)
AIC(m1, m2, m3, m4, m5)
summary(m1)

r.squaredGLMM(m4)
plot(gdDat$r2All, gdDat$betaG)
plot(gdDat$ntraits, gdDat$betaG)
plot(log(gdDat$dmean), gdDat$betaG)

#### Comparing e and c (Appendix 2) ####
x11(height=6, width=9)
par(mfrow=c(2,3), mar=c(4,4,2,2), oma=c(0,2,0,0), xpd=F)
plot(gdDat$betaT, gdDat$betaT_cond, las=1,
     xlim=c(-1.5, 2), ylim=c(-1.5, 2), pch=16,
     main="Original traits",
     xlab="",
     ylab="")
mtext("Slope for evolvabilities", 1, line=3)
mtext("Slope for cond. evolvabilities", 2, line=3)

lines(-10:10, -10:10)

plot(gdDat$betaD, gdDat$betaD_cond, las=1,
     xlim=c(0,4), ylim=c(0,4), pch=16,
     main="D-directions",
     xlab="",
     ylab="")
mtext("Slope for evolvabilities", 1, line=3)

lines(-10:10, -10:10)

plot(gdDat$betaP, gdDat$betaP_cond, las=1,
     xlim=c(-.5,3), ylim=c(-.5,3), pch=16,
     main="P-directions",
     xlab="",
     ylab="")
mtext("Slope for evolvabilities", 1, line=3)
lines(-10:10, -10:10)

plot(gdDat$r2T, gdDat$r2T_cond, las=1,
     xlim=c(0, 1), ylim=c(0, 1), pch=16,
     main="",
     xlab="",
     ylab="")
mtext(expression(paste(r^2, " for evolvabilities")), 1, line=3)
mtext(expression(paste(r^2, " for cond. evolvabilities")), 2, line=2.5)

lines(-10:10, -10:10)

plot(gdDat$r2D, gdDat$r2D_cond, las=1,
     xlim=c(0,1), ylim=c(0,1), pch=16,
     main="",
     xlab="",
     ylab="")
mtext(expression(paste(r^2, " for evolvabilities")), 1, line=3)

lines(-10:10, -10:10)

plot(gdDat$r2P, gdDat$r2P_cond, las=1,
     xlim=c(0,1), ylim=c(0,1), pch=16,
     main="",
     xlab="",
     ylab="")
mtext(expression(paste(r^2, " for evolvabilities")), 1, line=3)

lines(-10:10, -10:10)

