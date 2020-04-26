############################################################
#### Analyses of evolvability along divergence vectors #####
############################################################

rm(list=ls())

add2deltaList=function(pop=POPBASE[[s]]$popmeans[,1]){
  data.frame(species = out$species, g = out$gmat, ntraits = ncol(out$G),
             d=out$dmat, pop=pop, 
             emean=outdat$emean,
             emin=outdat$emin,
             emax=outdat$emax,
             cmean=outdat$cmean,
             div=outdat$div, edelta=outdat$edelta, cdelta=outdat$cdelta,
             theta=outdat$theta, row.names=NULL)
}

library(reshape2)
library(MCMCglmm)
library(evolvability)
source("code/prepareGD.R")
source("code/computeDelta.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")

# Species present in both databases
gsp = unlist(lapply(EVOBASE, function(x) x$Species))
dsp = unlist(lapply(POPBASE, function(x) x$Species))
names(dsp)
both_sp = unique(gsp[which(gsp %in% dsp)])
both_sp

deltaList = list()

#### Lobelia ####

# CERA G matrix, D = Caruso 2003 field
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means
colnames(means)==names(z0)

#outdat = computeDelta(G=out$G, means=means, z0=z0)
#evolvabilityMeans(out$G)
outdat = computeDelta2(G=out$G, means=means, z0=z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

# CERA G matrix, D=Caruso 2012 
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Means
s=out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1,3:6)] #Caruso 2012
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

# Krumm G matrix, D=Caruso 2003
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 3)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means #Caruso 2003
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

# Krumm G matrix, D=Caruso 2012 
out = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1,3:6)] #Caruso 2012
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList()

#### Aquilegia ####

# G = QFP1, D = Herlihy and Eckert 2007
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(means[2,])
means = means[-2,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-2,1])

# G = QFP1, D = Bartkowska et al. 2018
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# G = QLL3, D = Herlihy and Eckert 2007
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(means[3,]) #Second G mat, first D mat
means = means[-3,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-3,1])

# G = QLL3, D = Bartkowska et al. 2018
out = prepareGD(species = "Aquilegia_canadensis", gmatrix = 2, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means #Second D mat
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

#### Ipomopsis ####
# D = Caruso 2000
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 3)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0=EVOBASE[[out$gmat]]$Means[c(2:5)]
colnames(means)==names(z0)

outdat=computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# D = Caruso 2001
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 4)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(2,3,5)]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# D = Campbell 2018
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1:3)]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# D = Campbell 2019 unpubl.
out = prepareGD(species = "Ipomopsis_aggregata", gmatrix = 1, dmatrix = 2)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means[c(1:3, 5)]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

#### Brassica ####
gg=1
for(gg in 1:5){
  out = prepareGD(species = "Brassica_cretica", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s = out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(means[gg,])
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(G=out$G, means=means, z0=z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Turnera ####
out = prepareGD(species = "Turnera_ulmifolia", gmatrix = 1, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = EVOBASE[[out$gmat]]$Means
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

#### Clarkia ####
out = prepareGD(species = "Clarkia_dudleyana", gmatrix = 2, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(c(means[4,])) #Second G-matrix
means = means[-4,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-4,1])

#### Spergularia ####
for(gg in 1:4){
  out = prepareGD(species = "Spergularia_marina", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s=out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Solanum ####
for(gg in 1:3){
  out = prepareGD(species = "Solanum_carolinense", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s=out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  
  colnames(means)==names(z0)
  
  outdat=computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Fragaria ####

# G = Heterophrodite
out = prepareGD(species = "Fragaria_virginiana", gmatrix = 3, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = apply(means, 2, mean)
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

# G = Female
out = prepareGD(species = "Fragaria_virginiana", gmatrix = 4, dmatrix = 1)
out$gmat
out$dmat

# Means
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = apply(means, 2, mean)
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[,1])

##### Mimulus guttatus ####
gmats = c(5,7)

for(gg in 1:2){
out = prepareGD(species = "Mimulus_guttatus", gmatrix = gmats[gg], dmatrix = 1)
out$gmat
out$dmat

# Means
names(POPBASE)
s = out$dmat
means = POPBASE[[s]]$popmeans[,-1]
means = means[,match(colnames(out$D), names(means))]
z0 = unlist(c(means[gg,])) #First G-matrix
means = means[-gg,]
colnames(means)==names(z0)

outdat = computeDelta2(out$G, means, z0)

deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}


##### Mimulus micranthus ####
for(gg in 1:2){
  out = prepareGD(species = "Mimulus_micranthus", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Meangg
  s = out$dmat
  means = abs(POPBASE[[s]]$popmeans[,-1])
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

##### Nigella degenii ####
gg=1
for(gg in 1:2){
  out = prepareGD(species = "Nigella_degenii", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s = out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}

#### Lobelia ####
for(gg in 1:2){
  out = prepareGD(species = "Lobelia_siphilitica", gmatrix = gg, dmatrix = 1)
  out$gmat
  out$dmat
  
  # Means
  s = out$dmat
  means = POPBASE[[s]]$popmeans[,-1]
  means = means[,match(colnames(out$D), names(means))]
  z0 = unlist(c(means[gg,])) #First G-matrix
  means = means[-gg,]
  colnames(means)==names(z0)
  
  outdat = computeDelta2(out$G, means, z0)
  
  deltaList[[length(deltaList)+1]] = add2deltaList(pop=POPBASE[[s]]$popmeans[-gg,1])
}


#### Plotting deltaList ####

#save(deltaList, file="deltaList.RData")

load(file="deltaList.RData")

# Add data from orginal analyses
load(file="analyses/andersson_crepis/deltaDF.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/dalechampia/deltaDF_Tulum_CR.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/dalechampia/deltaDF_Tulum_MX.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/dalechampia/deltaDF_Tovar_MX.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/walter2018/deltaDF_Dune.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/walter2018/deltaDF_Head.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/walter2018/deltaDF_Table.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/walter2018/deltaDF_Wood.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/puentes2016/deltaDF_SPIT.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/puentes2016/deltaDF_STUC.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/puentes2016/deltaDF_STUS.RData")
deltaList[[length(deltaList)+1]]=deltaDF
load(file="analyses/puentes2016/deltaDF_VIS.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/colautti/deltaDF.RData")
deltaList[[length(deltaList)+1]]=deltaDF

load(file="analyses/delph_Silene/deltaDF.RData")
deltaList[[length(deltaList)+1]]=deltaDF

# Compile dataframe
library(plyr)
deltaDat = rbind.fill(deltaList)
dim(deltaDat)
tail(deltaDat)

for(i in 1: nrow(deltaDat)){
  if(deltaDat$theta[i]>90){
    deltaDat$theta[i]=180-deltaDat$theta[i]
  }
}

save(deltaDat, file="deltaDat.RData")

x11(height=4, width=10)
par(mfrow=c(1,3))
plot(deltaDat$edelta, deltaDat$theta, ylim=c(0, 90), las=1, pch=16, col="black",
     ylab="Angle between divergence vector and gmax",
     xlab="Evolvability along divergence vector")

plot(deltaDat$edelta/deltaDat$emean, deltaDat$theta, ylim=c(0,90), las=1, pch=16, col="black",
     ylab="Angle between divergence vector and gmax",
     xlab="Evolvability_delta/e_mean)")
abline(v=1, lty=2)

plot(deltaDat$edelta/deltaDat$emax, deltaDat$theta, ylim=c(0,90), las=1, pch=16, col="black",
     ylab="Angle between divergence vector and gmax",
     xlab="Evolvability_delta/e_max")

x11(height=4, width=10)
par(mfrow=c(1,3))
plot(deltaDat$edelta, deltaDat$theta, ylim=c(0, 90), las=1, pch=16, col="black",
     ylab="Angle between divergence vector and gmax",
     xlab="Evolvability along divergence vector")

plot(log(deltaDat$edelta/deltaDat$emean), deltaDat$theta, ylim=c(0,90), las=1, pch=16, col="black",
     ylab="Angle between divergence vector and gmax",
     xlab="log(e_delta/e_mean)")
abline(v=0, lty=2)

plot(deltaDat$edelta/deltaDat$emax, deltaDat$theta, ylim=c(0,90), las=1, pch=16, col="black",
     ylab="Angle between divergence vector and gmax",
     xlab="Evolvability_delta/e_max")


deltaDat[which(deltaDat$theta > 75 & (deltaDat$edelta/deltaDat$emax)>0.4),]

meanDat = ddply(deltaDat, .(species), summarize,
              emean = median(emean),
              emax = median(emax),
              emin = median(emin),
              cmean = median(cmean),
              edelta= median(edelta),
              cdelta = median(cdelta),
              div = median(div))

deltaDat=deltaDat[-which(deltaDat$sp=="Senecio_pinnatifolius"),]

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

#### Scaling between mean and max ####
plot(log10(deltaDat$emean), (deltaDat$edelta-deltaDat$emean)/(deltaDat$emax-deltaDat$emean), pch=16, las=1, main="Evolvability",
     col="lightgrey",
     xlab="",
     ylab="",
     xlim=c(-1, 1.5), ylim=c(-1,1))
#mtext(expression(paste("log[", e(Delta), "/", bar(e), "]")), 2, line=2)
mtext("Mean evolvability (%)", 1, line=2.5)
points(log10(meanDat$emean), (meanDat$edelta-meanDat$emean)/(meanDat$emax-meanDat$emean), pch=16)
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

deltaDat2=deltaDat
#deltaDat2 = deltaDat[-which(deltaDat$sp=="Mimulus micranthus"),]
#deltaDat2$cmean[which(deltaDat2$cmean<0.01)]=0.01
meanDat2=meanDat
#meanDat2 = meanDat2[-which(meanDat$sp=="Mimulus micranthus"),]
#deltaDat2$cdelta[which(deltaDat2$cdelta<0.01)]=0.01

plot(log10(deltaDat2$cmean), log(deltaDat2$cdelta/deltaDat2$cmean), pch=16, las=1, main="Conditional evolvability",
     col="lightgrey",
     xlab="",
     ylab="",
     xaxt="n",
     xlim=c(-2, 1), ylim=c(-2, 4))
axis(1, at=c(-2:1), labels = c(10^(-2:1)))
mtext(expression(paste("log[", c(Delta), "/", bar(c), "]")), 2, line=2)
mtext("Mean evolvability (%)", 1, line=2.5)
points(log10(meanDat2$cmean), log(meanDat2$cdelta/meanDat2$cmean), pch=16)
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

points(log10(median(meanDat2$cmean)), median(log(meanDat2$cdelta/meanDat2$cmean)), pch=16, col="blue")

plot(deltaDat$div, (deltaDat$cdelta-deltaDat$cmean)/(deltaDat$emax-deltaDat$cmean), pch=16, las=1, main="Evolvability",
     col="lightgrey",
     xlab="",
     ylab="",
     xlim=c(0, 120), ylim=c(-1,1))
#mtext(expression(paste("log[", e(Delta), "/", bar(e), "]")), 2, line=2)
mtext("Mean evolvability (%)", 1, line=2.5)
points(meanDat$div, (meanDat$cdelta-meanDat$cmean)/(meanDat$emax-meanDat$cmean), pch=16)

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

plot(deltaDat2$div, log(deltaDat2$cdelta/deltaDat2$cmean), pch=16, las=1,
     col="lightgrey",
     xlab="",
     ylab="",
     xlim=c(-5, 120), ylim=c(-2, 4))
points(meanDat2$div, log(meanDat2$cdelta/meanDat2$cmean), pch=16)
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

points(median(meanDat2$div), median(log(meanDat2$cdelta/meanDat2$cmean)), pch=16, col="blue")




newe = (deltaDat$edelta-deltaDat$emean)/(deltaDat$emax-deltaDat$emean)
for(i in 1:nrow(deltaDat)){
  if(newe[i]<0){
    newe[i] = ((deltaDat$edelta[i]-deltaDat$emin[i])/(deltaDat$emean[i]-deltaDat$emin[i]))-1
  }
}

newmeans = (meanDat$edelta-meanDat$emean)/(meanDat$emax-meanDat$emean)
for(i in 1:nrow(meanDat)){
  if(newe[i]<0){
    newe[i]=((meanDat$edelta[i]-meanDat$emin[i])/(meanDat$emean[i]-meanDat$emin[i]))-1
  }
}

par(mfrow=c(1,1))
plot(deltaDat$div, newe, pch=16, col="grey", ylim=c(-1,1), yaxt="n", ylab="", xlab="")
points(meanDat$div, newmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(e)), expression(e[max])), las=1)
mtext("Divergence from focal population (x100)", 1, line=2.5)

newe = (deltaDat$cdelta-deltaDat$cmean)/(deltaDat$emax-deltaDat$cmean)
for(i in 1:nrow(deltaDat)){
  if(newe[i]<0){
    newe[i] = ((deltaDat$cdelta[i]-emin[i])/(deltaDat$cmean[i]-emin[i]))-1
  }
}

newmeans = (meanDat$cdelta-meanDat$cmean)/(meanDat$emax-meanDat$cmean)
for(i in 1:nrow(meanDat)){
  if(newe[i]<0){
    newe[i]=((meanDat$cdelta[i]-meanDat$emin[i])/(meanDat$cmean[i]-meanDat$emin[i]))-1
  }
}

par(mfrow=c(1,1))
plot(deltaDat$div, newe, pch=16, col="grey", ylim=c(-1,1), yaxt="n", ylab="", xlab="")
points(meanDat$div, newmeans, pch=16)
abline(h=0)
axis(2, at=c(-1,0,1), labels=c(expression(e[min]), expression(bar(c)), expression(e[max])), las=1)
mtext("Divergence from focal population (x100)", 1, line=2.5)



#### Plotting individual studies from deltaList ####
par(mfrow=c(1,1))
plotDelta(44, lab.offset=0.05)

pdf("figs/delta_plots.pdf", family="Times", height=5, width=5)
for(i in 1:length(deltaList)){
  plotDelta(i, lab.offset=0.05)
}
dev.off()

plotDelta = function(index, lab.offset=0.05){
  plotData = deltaList[[index]]
  
  ylim = c(plotData$emin[1], plotData$emax[1])
  ylim = c(ylim[1]-abs(ylim[1]*.1), ylim[2]+abs(ylim[2]*.3))
  
  par(mar=c(4,4,5,4))
  plot(plotData$div, plotData$edelta, pch=16, ylim=ylim, las=1,
       xlab="", ylab="",
       main=paste(plotData$g[1],"/",
                  plotData$d[1]))
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
}


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


