##############################################
#### Lobelia plots for figure 3 ####
##############################################

rm(list=ls())

library(evolvability)
source("code/prepareGD.R")

load("data/EVOBASE.RData")
load("data/POPBASE.RData")
load("data/PMATBASE.RData")

# G = CERA, D = Caruso 2012
out = prepareGD(species="Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)

# Mean G-matrix
glist = list()
glist[[1]] = prepareGD(species = "Lobelia_siphilitica", gmatrix = 1, dmatrix = 2)$G
glist[[2]] = prepareGD(species = "Lobelia_siphilitica", gmatrix = 2, dmatrix = 2)$G
gmat = apply(simplify2array(glist), 1:2, mean)

# Load the error-corrected D matrix
load(file="analyses/adj_Dmats/Lobelia_Caruso2012.RData")
modD = matrix(apply(modDpost, 2, median), nrow=sqrt(ncol(modDpost)))

out$D = modD*100
out$G = gmat*100 #Mean G

# P matrix
plist = list()
ma = match(colnames(out$G), names(PMATBASE[["Lobelia siphilitica: CERA"]]$Means))

plist[[1]] = meanStdG(PMATBASE[["Lobelia siphilitica: CERA"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: CERA"]]$Means[ma])
plist[[2]] = meanStdG(PMATBASE[["Lobelia siphilitica: Krumm"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Krumm"]]$Means[ma])
plist[[3]] = meanStdG(PMATBASE[["Lobelia siphilitica: Reichelt"]]$P[ma, ma], 
                      PMATBASE[["Lobelia siphilitica: Reichelt"]]$Means[ma])
MeanP = apply(simplify2array(plist), 1:2, mean)

# Manual plot
dmat = out$D
gmat = out$G

g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e

p_ev = eigen(MeanP)$vectors
var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e

# Compute summary stats
mt = lm(log(diag(dmat))~log(diag(gmat)))
beta_t = summary(mt)$coef[2,1]
r2_t = summary(mt)$r.squared

mg = lm(log(var_d_g)~log(var_g_g))
beta_g = summary(mg)$coef[2,1]
r2_g = summary(mg)$r.squared

md = lm(log(var_d_d)~log(var_g_d))
beta_d = summary(md)$coef[2,1]
r2_d = summary(md)$r.squared

mp = lm(log(var_d_p)~log(var_g_p))
beta_p = summary(mp)$coef[2,1]
r2_p = summary(mp)$r.squared

# Plot
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(diag(gmat)), log10(diag(dmat)), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="", 
     ylab="",
     yaxt="n",
     xaxt="n",
     main="", las=1, pch=16, col="blue3")
mtext("Evolvability (%)", 1, line=2.5, cex=1)
mtext(expression(paste(italic(Lobelia), " ", italic(siphilitica))), line=0.5, cex=0.9)

points(log10(var_g_g), log10(var_d_g), pch=16)
points(log10(var_g_d), log10(var_d_d), pch=1)
points(log10(var_g_p), log10(var_d_p), pch=16, col="firebrick")

mean1 = mean(log10(c(diag(dmat), var_d_g, var_d_d)))
mean2 = mean(log10(c(diag(gmat), var_g_g, var_g_d)))
segments(x0=mean2-10, y0=mean1-10, x1=mean2+10, y1=mean1+10)

legend("bottomright", legend=c(paste0("Original traits (", round(100*r2_t, 1),"%)"),
                               paste0("G-directions (", round(100*r2_g, 1),"%)"),
                               paste0("D-directions (", round(100*r2_d, 1),"%)"),
                               paste0("P-directions (", round(100*r2_p, 1),"%)")),
       pch=c(16, 16, 1, 16), col=c("blue3", "black", "black", "firebrick"))

axis(1, at=seq(-.4, 1, .2), signif(10^seq(-.4, 1, .2), 2))

x3at = seq(-2, .5, .5)
x3 = exp(sqrt(((10^x3at)/100)*(2/pi)))
axis(2, at=x3at, signif(x3, 3), las=1)

dev.off() # For complex figure compiled from multiple R files

