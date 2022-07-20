################################################################################
##### Combine the results of the eigenvector and divergence vector analyses ####
################################################################################

rm(list=ls())
library(plyr)
library(lme4)

load(file="gdDatSE.RData")
load(file="deltaDat.RData")
load("data/EVOBASE.RData")
load("data/POPBASE.RData")

# Summary tables
names(gdDat)
GDTable = gdDat
GDTable[,c(9:13, 15:34)] = signif(GDTable[,c(9:13, 15:34)], 3)
#write.csv2(GDTable, file = "GDTable.csv", dec=".")

names(deltaDat)
DeltaTable = deltaDat
#write.csv2(DeltaTable, file = "DeltaTable.csv", dec=".")

# Combine dataframes
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
                        betaP = median(betaP, na.rm=T),
                        betaP_cond = median(betaP_cond, na.rm=T),
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
               betaD_SE = betaD_SE,
               r2D = median(r2D),
               betaP = median(betaP, na.rm=T),
               betaP_SE = betaP_SE,
               betaP_cond = median (betaP_cond, na.rm=T),
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

# Start of analyses ####
table(comb2$traitgroups)
table(comb2$dims)
table(comb2$ndims)

# Figures for Appendix 1 (multivariate scaling relationships) ####

# Trait groups
x11(height=5, width=8)
par(mfrow=c(1,2))
par(mar=c(8,4,2,2))

medians = tapply(comb2$betaG, comb2$traitgroups, median, na.rm=T)
comb2$traitgroups = factor(comb2$traitgroups, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$traitgroups), comb2$betaG, xaxt="n", xlab="", las=1,
     ylab="Slope for G-directions")
axis(1, at=1:4, labels=rep("", 4))
text(1:4, par("usr")[3] - .1, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$traitgroups), xpd = TRUE)

medians = tapply(comb2$edelta/comb2$emean.y, comb2$traitgroups, median, na.rm=T)
comb2$traitgroups = factor(comb2$traitgroups, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$traitgroups), comb2$edelta/comb2$emean.y, xaxt="n", xlab="", las=1,
     ylab="")
mtext(expression(paste("Prop evol along div vector  "," [", e(Delta), "/", bar(e), "]")), 2, line=2.5)

axis(1, at=1:4, labels=rep("", 4))
text(1:4, par("usr")[3] - .25, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$traitgroups), xpd = TRUE)

# Mating systems

# Get mating systems from EVOBASE
ms = NULL
for(i in 1:nrow(comb2)){
ms[i] = paste0(EVOBASE[[as.character(comb2$g[i])]]$MS, collapse = "+")
}
ms = factor(ms, levels=c("S", "M", "O"))
comb2$ms = ms
table(comb2$ms)

# Plot
x11(height=5, width=8)
par(mfrow=c(1,2))
plot(ms, comb2$betaG, las=1, xlab="", ylab="Slope for G-directions")
plot(ms, comb2$edelta/comb2$emean.y, las=1, xlab="", ylab="")
mtext(expression(paste("Prop evol along div vector   "," [", e(Delta), "/", bar(e), "]")), 2, line=2.5)

# Study environments

# Get study environments (for the divergence data) from POPBASE
env = NULL
for(i in 1:nrow(comb2)){
  env[[i]] = paste0(POPBASE[[as.character(comb2$d[i])]]$Env, collapse="+")
}
env = factor(env, levels=c("field", "common_garden", "greenhouse"))
comb2$env = env
table(comb2$env)

# Plot
x11(height=5*1.3, width=8*1.3)
par(mfrow=c(1,2))
plot(env, comb2$betaG, las=1, xlab="", ylab="Slope for G-directions")
plot(env, comb2$edelta/comb2$emean.y, las=1, xlab="", ylab="")
mtext(expression(paste("Prop evol along div vector   "," [", e(Delta), "/", bar(e), "]")), 2, line=2.5)

#### Figures for Appendix 3 (comparing GD and delta analyses) ####
x11(height=5, width=9)
par(mfrow=c(1,2))
plot(log(comb2$edelta/comb2$emean.x), comb2$betaG, pch=16, col="lightgrey", las=1,
     xlab="", ylab="Slope of G-directions")
points(log(comb$edelta/comb$emean.x), comb$betaG, pch=16)
mtext(expression(paste("Prop evol along div vector "," (log "," [", e(Delta), "/", bar(e), "])")), 1, line=2.5)
abline(v=0, lty=2)
abline(h=1, lty=2)

plot(log(comb2$cdelta/comb2$cmean.x), comb2$betaG, pch=16, col="lightgrey", las=1,
     xlab="", ylab="Slope of G-directions")
points(log(comb$cdelta/comb$cmean.x), comb$betaG, pch=16)
mtext(expression(paste("Prop cond evol along div vector "," (log "," [", c(Delta), "/", bar(c), "])")), 1, line=2.5)
abline(v=0, lty=2)
abline(h=1, lty=2)

# Angles etc.
comb2$r2M = apply(subset(comb2, select =c("r2G", "r2D")), 1, mean, na.rm=T)

x11(height=11, width=12)
par(mfrow=c(2,2), mar = c(4,4,1,1))
plot(comb2$r2All, comb2$thetaGD, pch=16, col="lightgrey", las=1, xlim=c(0,1), ylim=c(0,90),
     xlab="", ylab="")
points(comb$r2All, comb$thetaGD, pch=16)
mtext(expression(paste("Overall ", r^2, " of scaling relationship")), 1, line=2.5, cex=.9)
mtext(expression(paste("Angle between  ", g[max], " and ", d[max])), 2, line=2.5, cex=.9)

plot(comb2$thetaGD, comb2$betaG, pch=16, col="lightgrey", las=1,
     ylab="", xlab="")
points(comb$thetaGD, comb$betaG, pch=16)
mtext(expression(paste("Angle between ", g[max], " and ", d[max])), 1, line=2.5, cex=.9)
mtext("Slope for G-directions", 2, line=2.5, cex=.9)

plot(log(comb2$edelta/comb2$emean.x), comb2$thetaGD, pch=16, col="lightgrey", las=1,
     xlab="", ylab="")
points(log(comb$edelta/comb$emean.x), comb$thetaGD, pch=16)
mtext(expression(paste("Prop evol along div vector ","(log"," [", e(Delta), "/", bar(e), "])")), 1, line=2.5, cex=.9)
mtext(expression(paste("Angle between ", g[max], " and ", d[max])), 2, line=2.5, cex=.9)
abline(v=0, lty=2)

plot(comb2$edelta/comb2$emax, comb2$thetaGD, pch=16, col="lightgrey", las=1,
     xlab="", ylab="")
points(comb$edelta/comb$emax, comb$thetaGD, pch=16)
mtext(expression(paste("Angle between ", g[max], " and ", d[max])), 2, line=2.5, cex=.9)
mtext(expression(paste("Prop of max evol along div vector ","[", e(Delta), "/", e[max], "]")), 1, line=2.5, cex=.85)

# Figures for Appendix 4 (trait dimensions) ####
x11(height=11, width=12)
par(mfrow=c(2,2), mar=c(3,4,4,2))

plot(factor(comb2$ndims), comb2$betaG, xlab="", ylab="", las=1)
mtext("Number of trait dimensions", 1, line=2.5, cex=.9)
mtext("Slope for G-directions", 2, line=2.5, cex=.9)

plot(factor(comb2$ndims), comb2$edelta/comb2$emean.y, xlab="", ylab="", las=1)
mtext("Number of trait dimensions", 1, line=2.5, cex=.9)
mtext(expression(paste("Prop evol along div vector  ","[", e(Delta), "/", bar(e), "]")), 2, line=2, cex=.9)

par(mar=c(8,4,2,2))
medians = tapply(comb2$betaG, comb2$dims, median, na.rm=T)
comb2$dims = factor(comb2$dims, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$dims), comb2$betaG, xaxt="n", xlab="", las=1, ylab="")
mtext("Slope for G-directions", 2, line=2.5, cex=.9)

axis(1, at=1:10, labels=rep("", 10))
text(1:10, par("usr")[3] - .15, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$dims), xpd = TRUE)

medians = tapply(comb2$edelta/comb2$emean.y, comb2$dims, median, na.rm=T)
comb2$dims = factor(comb2$dims, levels=names(sort(medians, decreasing=F)))
plot(factor(comb2$dims), comb2$edelta/comb2$emean.y, xaxt="n", xlab="", las=1,
     ylab="")
mtext(expression(paste("Prop evol along div vector  ", "[", e(Delta), "/", bar(e), "]")), 2, line=2, cex=.9)

axis(1, at=1:10, labels=rep("", 10))
text(1:10, par("usr")[3] - .4, srt = 45, adj = 1, cex=1,
     labels = levels(comb2$dims), xpd = TRUE)

#### END OF FINAL ANALYSES ####