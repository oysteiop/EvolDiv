################################################################################
##### Postprocessing and figures for divergence-vector analyses ####
################################################################################

rm(list=ls())

load(file="deltaListFull.RData")
load(file="deltaDat.RData")

# Misc plots for supporting information

# Plots of delta vs. angles ####
x11(height=3.5, width=10)
par(mfrow=c(1,3), mar=c(4,4,2,2))
plot(deltaDat$edelta, deltaDat$theta, ylim=c(0, 90), las=1, pch=16, col="black",
     ylab="",
     xlab="")

mtext(expression(paste("Angle btw div vector and  ",g[max])), 2, line=2.5, cex=.9)
mtext(expression(paste("Evolvability along div vector [", e(Delta),"]")), 1, line=2.5, cex=.9)

plot(deltaDat$edelta/deltaDat$emean, deltaDat$theta, ylim=c(0,90), las=1, pch=16, col="black",
     ylab="",
     xlab="")
abline(v=1, lty=2)
#mtext(expression(paste("Angle btw div vector and  ",g[max])), 2, line=2.5)
mtext(expression(paste("Prop evol along div vector", " [", e(Delta), "/", bar(e), "]")), 1, line=2.5, cex=.9)

plot(deltaDat$edelta/deltaDat$emax, deltaDat$theta, ylim=c(0,90), las=1, pch=16, col="black",
     ylab="",
     xlab="")
mtext(expression(paste("Prop of max evol along div vector "," [", e(Delta), "/", e[max], "]")), 1, line=2.5, cex=.9)


deltaDat[which(deltaDat$theta > 75 & (deltaDat$edelta/deltaDat$emax)>0.4),]

meanDat = ddply(deltaDat, .(species), summarize,
                emean = median(emean),
                emax = median(emax),
                emin = median(emin),
                cmean = median(cmean),
                edelta= median(edelta),
                cdelta = median(cdelta),
                div = median(div))

#deltaDat=deltaDat[-which(deltaDat$sp=="Senecio_pinnatifolius"),]

# Summary figure for divergence vector analyses ####
x11(height=6.5, width=8)
par(mfrow=c(2,2), mar=c(4,4,2,4))

plot(log10(deltaDat$emean), log(deltaDat$edelta/deltaDat$emean), pch=16, las=1, main="Evolvability",
     col="lightgrey",
     xlab="",
     ylab="",
     xaxt="n",
     xlim=c(-1, 1.5), ylim=c(-2, 2))
axis(1, at=c(-1:1), labels = c(10^(-1:1)))
mtext(expression(paste("log[", e(Delta), "/", bar(e), "]")), 2, line=2)
mtext("Mean evolvability (%)", 1, line=2.5)
points(log10(meanDat$emean), log(meanDat$edelta/meanDat$emean), pch=16)
abline(h=0, lty=1)
abline(h=log(1/2), lty=2)
abline(h=log(1/3), lty=2)
abline(h=log(2), lty=2)
abline(h=log(3), lty=2)
x=1.7
text(x=x, y=0, labels=expression(paste(e(Delta),"=",bar(e))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(2), labels=expression(paste(e(Delta),"=2",bar(e))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(3), labels=expression(paste(e(Delta),"=3",bar(e))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(1/2), labels=expression(paste(e(Delta),"=",frac(1,2),bar(e))), xpd=T, cex=.8, adj=0)

points(log10(median(meanDat$emean)), median(log(meanDat$edelta/meanDat$emean)), pch=16, col="blue")

plot(log10(deltaDat$cmean), log(deltaDat$cdelta/deltaDat$cmean), pch=16, las=1, main="Conditional evolvability",
     col="lightgrey",
     xlab="",
     ylab="",
     xaxt="n",
     xlim=c(-2, 1), ylim=c(-2, 4))
axis(1, at=c(-2:1), labels = c(10^(-2:1)))
mtext(expression(paste("log[", c(Delta), "/", bar(c), "]")), 2, line=2)
mtext("Mean evolvability (%)", 1, line=2.5)
points(log10(meanDat$cmean), log(meanDat$cdelta/meanDat$cmean), pch=16)
abline(h=0, lty=1)
abline(h=log(1/2), lty=2)
abline(h=log(1/3), lty=2)
abline(h=log(2), lty=2)
abline(h=log(3), lty=2)
x=1.2
text(x=x, y=0, labels=expression(paste(c(Delta),"=",bar(c))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(2), labels=expression(paste(c(Delta),"=2",bar(c))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(3), labels=expression(paste(c(Delta),"=3",bar(c))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(1/2), labels=expression(paste(c(Delta),"=",frac(1,2),bar(c))), xpd=T, cex=.8, adj=0)

points(log10(median(meanDat$cmean)), median(log(meanDat$cdelta/meanDat$cmean)), pch=16, col="blue")

plot(deltaDat$div, log(deltaDat$edelta/deltaDat$emean), pch=16, las=1,
     col="lightgrey",
     xlab="",
     ylab="",
     xlim=c(-5, 120), ylim=c(-2, 2))
points(meanDat$div, log(meanDat$edelta/meanDat$emean), pch=16)
mtext(expression(paste("log[", e(Delta), "/", bar(e), "]")), 2, line=2)
mtext("Divergence from focal population (x100)", 1, line=2.5)
abline(h=0, lty=1)
abline(h=log(1/2), lty=2)
abline(h=log(1/3), lty=2)
abline(h=log(2), lty=2)
abline(h=log(3), lty=2)
x=127
text(x=x, y=0, labels=expression(paste(e(Delta),"=",bar(e))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(2), labels=expression(paste(e(Delta),"=2",bar(e))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(3), labels=expression(paste(e(Delta),"=3",bar(e))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(1/2), labels=expression(paste(e(Delta),"=",frac(1,2),bar(e))), xpd=T, cex=.8, adj=0)

points(median(meanDat$div), median(log(meanDat$edelta/meanDat$emean)), pch=16, col="blue")

plot(deltaDat$div, log(deltaDat$cdelta/deltaDat$cmean), pch=16, las=1,
     col="lightgrey",
     xlab="",
     ylab="",
     xlim=c(-5, 120), ylim=c(-2, 4))
points(meanDat$div, log(meanDat$cdelta/meanDat$cmean), pch=16)
mtext(expression(paste("log[", c(Delta), "/", bar(c), "]")), 2, line=2)
mtext("Divergence from focal population (x100)", 1, line=2.5)
abline(h=0, lty=1)
abline(h=log(1/2), lty=2)
abline(h=log(1/3), lty=2)
abline(h=log(2), lty=2)
abline(h=log(3), lty=2)
x=127
text(x=x, y=0, labels=expression(paste(c(Delta),"=",bar(c))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(2), labels=expression(paste(c(Delta),"=2",bar(c))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(3), labels=expression(paste(c(Delta),"=3",bar(c))), xpd=T, cex=.8, adj=0)
text(x=x, y=log(1/2), labels=expression(paste(c(Delta),"=",frac(1,2),bar(c))), xpd=T, cex=.8, adj=0)

points(median(meanDat$div), median(log(meanDat$cdelta/meanDat$cmean)), pch=16, col="blue")

# Alternative summary figure with scaling between min, mean, and max ####

#deltaDat = deltaDat[deltaDat$ndims==1,]

x11(height=4, width=8)
par(mfrow=c(1,2), mar=c(4,4,2,2))

newe = (deltaDat$edelta-deltaDat$emean)/(deltaDat$emax-deltaDat$emean)
for(i in 1:nrow(deltaDat)){
  if(newe[i]<0){
    newe[i] = ((deltaDat$edelta[i]-deltaDat$emin[i])/(deltaDat$emean[i]-deltaDat$emin[i]))-1
  }
}

newmeans = (meanDat$edelta-meanDat$emean)/(meanDat$emax-meanDat$emean)
for(i in 1:nrow(meanDat)){
  if(newmeans[i]<0){
    newmeans[i]=((meanDat$edelta[i]-meanDat$emin[i])/(meanDat$emean[i]-meanDat$emin[i]))-1
  }
}

plot(deltaDat$div, newe, pch=16, col="black", ylim=c(-1,1), yaxt="n", ylab="", xlab="",
     main=" (e) All studies (evolvability)")
#points(meanDat$div, newmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(e)), expression(e[max])), las=1)
mtext("Divergence from focal population (x100)", 1, line=2.5)

newc = (deltaDat$cdelta-deltaDat$cmean)/(deltaDat$emax-deltaDat$cmean)
for(i in 1:nrow(deltaDat)){
  if(newc[i]<0){
    newc[i] = ((deltaDat$cdelta[i]-deltaDat$emin[i])/(deltaDat$cmean[i]-deltaDat$emin[i]))-1
  }
}

newcmeans = (meanDat$cdelta-meanDat$cmean)/(meanDat$emax-meanDat$cmean)
for(i in 1:nrow(meanDat)){
  if(newcmeans[i]<0){
    newcmeans[i]=((meanDat$cdelta[i]-meanDat$emin[i])/(meanDat$cmean[i]-meanDat$emin[i]))-1
  }
}

plot(deltaDat$div, newc, pch=16, col="grey", ylim=c(-1,1), yaxt="n", ylab="", xlab="",
     main=" (f) All studies (conditional evolvability)")
#points(meanDat$div, newcmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(c)), expression(e[max])), las=1)
mtext("Divergence from focal population (x100)", 1, line=2.5)

#### Function for plotting individual studies from deltaList ####
plotDelta = function(index, title=paste(plotData$g[1], "/", plotData$d[1]), 
                     xlab=T, ylab=T,
                     cex=.8, lab.offset=0.05){
  plotData = deltaList[[index]]
  
  ylim = c(plotData$emin[1], plotData$emax[1])
  ylim = c(ylim[1]-abs(ylim[1]*.1), ylim[2]+abs(ylim[2]*.3))
  
  #par(mar=c(4,4,5,4))
  plot(plotData$div, plotData$edelta, pch=16, ylim=ylim, las=1,
       xlab="", ylab="",
       main=title)
  if(xlab) {mtext("Divergence from focal population", 1, line=2.5, cex=.9)}
  if(ylab) {mtext("Evolvability along divergence vector (%)", 2, line=2.5, cex=.9)}
  points(plotData$div, plotData$cdelta, pch=16, col="grey")
  
  abline(h=plotData$emean, lty=2)
  abline(h=plotData$cmean, lty=2, col="grey")
  abline(h=plotData$emin, lty=1)
  abline(h=plotData$emax, lty=1)
  abline(h=plotData$edrift, lty=2, col="red")
  x = max(plotData$div)+(max(plotData$div)-min(plotData$div))*lab.offset
  text(x=x, y=plotData$emean[1], labels=expression(bar(e)), xpd=T, cex=cex, adj=0)
  text(x=x, y=plotData$emin[1], labels=expression(e[min]), xpd=T, cex=cex, adj=0)
  text(x=x, y=plotData$emax[1], labels=expression(e[max]), xpd=T, cex=cex, adj=0)
  text(x=x, y=plotData$cmean[1], labels=expression(bar(c)), xpd=T, cex=cex, adj=0)
  
  legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")
}

unlist(lapply(deltaList, function(x) unique(x$species)))

# Figure 4 ####
cairo_pdf("pubfigs/edeltaFig.pdf", height=6, width=10, fam="Times")
#x11(height=6, width=10)
par(mfrow=c(2,3), mar=c(4, 4, 2, 3))

plotDelta(37, lab.offset=0.07, cex=1, xlab=F, title = "")
mtext(expression(paste("", italic(Crepis)," ", italic(tectorum))), 3, line=0.5, cex=.8)
plotDelta(38, lab.offset=0.07, cex=1, xlab=F, ylab=F, title = "")
mtext(expression(paste("", italic(Dalechampia)," ", italic(scandens))), 3, line=0.5, cex=.8)

plot(deltaDat$div, newe, pch=16, col="black", ylim=c(-1,1), yaxt="n", ylab="", xlab="",
     main="")
mtext(expression("All studies (evolvability)"), 3, line=0.5, cex=.9)

#points(meanDat$div, newmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(e)), expression(e[max])), las=1)
#mtext("Divergence from focal population (x100)", 1, line=2.5)

plotDelta(45, lab.offset=0.07, cex=1, title = "")
mtext(expression(paste("", italic(Arabidopsis)," ", italic(lyrata))), 3, line=0.5, cex=.8)

plotDelta(2, lab.offset=0.07, cex=1, ylab=F, title = "")
mtext(expression(paste("", italic(Lobelia)," ", italic(siphilitica))), 3, line=0.5, cex=.8)

plot(deltaDat$div, newc, pch=16, col="grey", ylim=c(-1,1), yaxt="n", ylab="", xlab="",
     main="")
mtext(expression("All studies (cond. evolvability)"), 3, line=0.5, cex=.9)

#points(meanDat$div, newcmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(c)), expression(e[max])), las=1)
mtext("Divergence from focal population", 1, line=2.5, cex=.9)

dev.off()




pdf("figs/delta_plots.pdf", family="Times", height=5, width=5)
for(i in 1:length(deltaList)){
  plotDelta(i, lab.offset=0.05)
}
dev.off()



#### Plotting individual systems from deltaList ####
par(mfrow=c(1,1))

lab.offset=0.05

unique(deltaDat$g)

#D. scandensoides
plotData = deltaDat[deltaDat$g=="Tulum",]

plot(plotData$edelta, plotData$theta)


ylim = c(plotData$emin[1], plotData$emax[1])
ylim = c(ylim[1]-abs(ylim[1]*.1), ylim[2]+abs(ylim[2]*.3))

par(mar=c(4,4,5,4))
plot(plotData$div, plotData$edelta, pch=16, ylim=ylim, las=1, col=c("black", "blue")[as.numeric(factor(plotData$d))],
     xlab="", ylab="",
     main="D. scandensoides MX + CR")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(plotData$div, plotData$cdelta, pch=16, col="grey")

abline(h=plotData$emean, lty=2)
abline(h=plotData$cmean, lty=2, col="grey")
abline(h=plotData$emin, lty=1)
abline(h=plotData$emax, lty=1)
x = max(plotData$div)+(max(plotData$div)-min(plotData$div))*lab.offset
text(x=x, y=plotData$emean[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emin[1], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emax[1], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$cmean[1], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)

legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

#Lobelia
plotData = deltaDat[deltaDat$g=="Lobelia siphilitica: CERA",]
plotData$traits

ylim = c(plotData$emin[1], plotData$emax[1])
ylim = c(ylim[1]-abs(ylim[1]*.1), ylim[2]+abs(ylim[2]*.3))

par(mar=c(4,4,5,4))
plot(plotData$div, plotData$edelta, pch=16, ylim=ylim, las=1, col=c("black", "blue")[as.numeric(factor(plotData$d))],
     xlab="", ylab="",
     main="Lobelia siphilitica: CERA")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(plotData$div, plotData$cdelta, pch=16, col="grey")

abline(h=plotData$emean, lty=2)
abline(h=plotData$cmean, lty=2, col="grey")
abline(h=plotData$emin, lty=1)
abline(h=plotData$emax, lty=1)
x = max(plotData$div)+(max(plotData$div)-min(plotData$div))*lab.offset
text(x=x, y=plotData$emean[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emin[1], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emax[1], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$cmean[1], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)

legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

#Ipomopsis
plotData = deltaDat[deltaDat$g=="Ipomopsis aggregata: Vera Falls",]
plotData$traits

ylim = c(plotData$emin[1], plotData$emax[1])
ylim = c(ylim[1]-abs(ylim[1]*.1), ylim[2]+abs(ylim[2]*.3))

par(mar=c(4,4,5,4))
plot(plotData$div, plotData$edelta, pch=16, ylim=ylim, las=1, col=c("black", "blue")[as.numeric(factor(plotData$d))],
     xlab="", ylab="",
     main="Lobelia siphilitica: CERA")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(plotData$div, plotData$cdelta, pch=16, col="grey")

abline(h=plotData$emean, lty=2)
abline(h=plotData$cmean, lty=2, col="grey")
abline(h=plotData$emin, lty=1)
abline(h=plotData$emax, lty=1)
x = max(plotData$div)+(max(plotData$div)-min(plotData$div))*lab.offset
text(x=x, y=plotData$emean[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emin[1], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emax[1], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$cmean[1], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)

legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")

#Aquilegia
plotData = deltaDat[deltaDat$g=="Aquilegia canadensis: QFP1",]
plotData$traits

ylim = c(plotData$emin[1], plotData$emax[1])
ylim = c(ylim[1]-abs(ylim[1]*.1), ylim[2]+abs(ylim[2]*.3))

par(mar=c(4,4,5,4))
plot(plotData$div, plotData$edelta, pch=16, ylim=ylim, las=1, col=c("black", "blue")[as.numeric(factor(plotData$d))],
     xlab="", ylab="",
     main="Aquilegia canadensis: QFP1")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(plotData$div, plotData$cdelta, pch=16, col="grey")

abline(h=plotData$emean, lty=2)
abline(h=plotData$cmean, lty=2, col="grey")
abline(h=plotData$emin, lty=1)
abline(h=plotData$emax, lty=1)
x = max(plotData$div)+(max(plotData$div)-min(plotData$div))*lab.offset
text(x=x, y=plotData$emean[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emin[1], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$emax[1], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=plotData$cmean[1], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)

legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")


# All the studies ####
nG = NULL
nD = NULL
for(s in c(1:length(both_sp))){
  species = both_sp[s]  
  nG[s] = length(EVOBASE[which(unlist(lapply(EVOBASE, function(x) x$Species))==species)])
  nD[s] = length(POPBASE[which(unlist(lapply(POPBASE, function(x) x$Species))==species)])
}
cbind(both_sp, nG, nD)

reslist = list()
for(s in 1:length(both_sp)){
  for(g in 1:nG[s]){
    for(d in 1:nD[s]){
      res = prepareGD(species = both_sp[s], gmatrix = g, dmatrix = d)[c(1:2, 5:18)]
      reslist[length(reslist)+1] = as.data.frame(unlist(res)[c(1:16)])
    }
  }
}
length(reslist)
reslist[[20]]

resmat = as.data.frame(matrix(NA, ncol=16, nrow=length(reslist)))
for(i in 1:length(reslist)){
  resmat[i,1:2]=as.character(reslist[[i]][1:2])
  resmat[i,3:13]=round(as.numeric(as.character(reslist[[i]][3:13])),2)
  resmat[i,14:16]=as.numeric(as.character(reslist[[i]][14:16]))
}

colnames(resmat)=c("D", "G", "npop", "theta", "betaG", "r2G", "betaD", "r2D","betaT","r2T","i_mean","e_mean","d_mean","nBetaG","nBetaD","nBetaT")

for(i in 1:nrow(resmat)){
  if(resmat$theta[i]>=90){
    resmat$theta[i] = 180 - resmat$theta[i]
  }
}

names(resmat)
#resmat = resmat[,c(2,1,16,3, 12,13,4,5,6)]
head(resmat)
View(resmat)

plot(resmat$theta, resmat$betaG, pch=16)
plot(resmat$theta, resmat$r2G, pch=16)

resmat = resmat[resmat$npop>2,]
resmat = resmat[resmat$nBetaG>2,]
resmat = resmat[resmat$nBetaD>2,]
resmat = resmat[resmat$nBetaT>2,]

mean(resmat$betaG)
hist(resmat$r2G)
hist(resmat$betaG)

plot(resmat$betaT, resmat$betaG)
plot(resmat$r2T, resmat$r2G)

plot(resmat$theta, resmat$betaG)

plot(resmat$npop, resmat$betaG)
plot(resmat$npop, resmat$r2G, cex=resmat$nBetaG*.5)

plot(resmat$i_mean, resmat$r2G, cex=resmat$nBetaG*.5, xlim=c(0,1))

#### Error-corrected D original code ####

drop = which(is.na(rowSums(means[-1])))
if(length(drop)>0){
  means=means[-drop]
  eV=ev[-drop]
}

#Set prior
n = ncol(means)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1))

means[,1:ncol(means)] = apply(means[,1:ncol(means)], 2, function(x) x*10)
mev = melt(eV)$value*100
data = means
vars = paste0(colnames(means), collapse=", ")
vars
vars = paste0("c(",vars,") ~-1+trait")

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

mod<-MCMCglmm(as.formula(noquote(vars)),
              rcov = ~us(trait):units,
              mev = mev,
              data = means, 
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)

modD = matrix(apply(mod$VCV, 2, median)[2:(1+(ncol(means))^2)], nrow=n)
colnames(modD) = rownames(modD) = colnames(means)
modD = meanStdG(modD, colMeans(means))
round(modD, 3)

index=1




##### Holcus lanatus ####

# G matrices very poorly estimated, excluded
names(EVOBASE)
m1 = EVOBASE[["Holcus lanatus: Improved"]]$Means[c(2,5,6,7)]
m2 = EVOBASE[["Holcus lanatus: Traditional"]]$Means[c(1,4,5,6)]
names(m1)==names(m2)

g1 = droptraits(EVOBASE[["Holcus lanatus: Improved"]]$G)
g1 = meanStdG(g1, m1)*100
c(names(m1)==colnames(g1))

delta = log(m2)-log(m1)
delta = delta/sqrt(sum(delta^2)) #Unit-length

evolvabilityMeans(g1)
cov2cor(g1)
evolvabilityBeta(g1, delta)$e
evolvabilityBeta(g1, delta)$c
evolvabilityBeta(g1, delta)$e/evolvabilityMeans(g1)[1]

m1 = EVOBASE[["Holcus lanatus: Improved"]]$Means[c(3,5,6,7)]
m2 = EVOBASE[["Holcus lanatus: Traditional"]]$Means[c(2,4,5,6)]
names(m1)==names(m2)

g2 = droptraits(EVOBASE[["Holcus lanatus: Traditional"]]$G)
g2 = meanStdG(g2, m2)*100
c(names(m1)==colnames(g2))

delta = log(m2)-log(m1)
delta = delta/sqrt(sum(delta^2)) #Unit-length

cov2cor(g2)
evolvabilityMeans(g2)
evolvabilityBeta(g2, delta)$e
evolvabilityBeta(g2, delta)$e/evolvabilityMeans(g2)[1]

# Fenster & Carr 1997
m1 = EVOBASE[["Mimulus guttatus: S II"]]$Means
m2 = EVOBASE[["Mimulus guttatus: T II"]]$Means
names(m1)==names(m2)

#g1 = droptraits(EVOBASE[["Mimulus guttatus: S II"]]$G) #Exclude because a value of 0 precludes calculations
g2 = droptraits(EVOBASE[["Mimulus guttatus: T II"]]$G)
#g1 = meanStdG(g1, m1)*100
g2 = meanStdG(g2, m2)*100
c(names(m1)==colnames(g1))

delta = log(m1)-log(m2)
div = mean(abs(delta))*100
delta = delta/sqrt(sum(delta^2)) #Unit-length

deltaList[[length(deltaList)+1]] = data.frame(sp="Mimulus guttatus", g="Mimulus guttatus: T II", traits=ncol(g2), 
                                              d=s, 
                                              pop="S", 
                                              emean=evolvabilityMeans(g2)[1],
                                              emin=evolvabilityMeans(g2)[2],
                                              emax=evolvabilityMeans(g2)[3],
                                              cmean=evolvabilityMeans(g2)[4],
                                              div=div, edelta=evolvabilityBeta(g2, delta)$e, 
                                              cdelta=evolvabilityBeta(g2, delta)$c,
                                              theta=acos(t(eigen(g2)$vectors[,1]) %*% delta)*(180/pi),
                                              row.names=NULL)

