################################################################################
##### Postprocessing and figures for divergence-vector analyses ####
################################################################################

rm(list=ls())

library(plyr)
library(lme4)
library(glmmTMB)

load(file="deltaListFullSE.RData")
load(file="deltaDatSE.RData")

# Meta-analysis of proportional evolvability along divergence vectors
eratio = deltaDat$edelta/deltaDat$emean
cratio = deltaDat$cdelta/deltaDat$cmean

sum(eratio>1)/sum(eratio>-Inf)*100
sum(cratio>1)/sum(cratio>-Inf)*100

hist(eratio)
hist(cratio)
hist(log(eratio))
hist(log(cratio))

mean(eratio)
mean(cratio)
median(eratio)
median(cratio)

deltaDat$pop = paste0(deltaDat$species, deltaDat$pop)

# Meta-analysis
ww = 1/(1+(deltaDat$edeltaSE^2)/(eratio^2))

m = glmmTMB(log(eratio) ~ 1 + (1|g) + (1|pop), weights = ww, REML = T, data=deltaDat)

summary(m)
exp(0.48328)

ww = 1/(1+(deltaDat$cdeltaSE^2)/(cratio^2))
m = glmmTMB(log(cratio) ~ 1 + (1|g) + (1|pop),  weights = ww, REML = T, data=deltaDat)

summary(m)
exp(0.8281)

# Misc plots for supporting information ####

# Plots of delta vs. angles
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

#deltaDat[which(deltaDat$theta > 75 & (deltaDat$edelta/deltaDat$emax)>0.4),]

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
     main="All studies (evolvability)")
#points(meanDat$div, newmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(e)), expression(e[max])), las=1)
mtext("Divergence from focal population", 1, line=2.5)

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
     main="All studies (conditional evolvability)")
#points(meanDat$div, newcmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(c)), expression(e[max])), las=1)
mtext("Divergence from focal population (x100)", 1, line=2.5)

#### Function for plotting individual studies from deltaList ####
plotDelta = function(index, title=paste(plotData$g[1], "/", plotData$d[1]), 
                     xlab=T, ylab=T, clab=T,
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
  if(clab){text(x=x, y=plotData$cmean[1], labels=expression(bar(c)), xpd=T, cex=cex, adj=0)}
  
  legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")
}

unlist(lapply(deltaList, function(x) unique(x$species)))

# Figure 4 ####
cairo_pdf("pubfigs/edeltaFig.pdf", height=6, width=10, fam="Times")
#x11(height=6, width=10)
par(mfrow=c(2,3), mar=c(4, 4, 2, 3))

plotDelta(37, lab.offset=0.07, cex=1, xlab=F, clab=F, title = "")
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
