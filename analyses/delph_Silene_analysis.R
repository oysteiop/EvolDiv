####################################################
#### Silene latifolia, data from Yu et al. 2011 ####
####################################################
rm(list=ls())
library(evolvability)
library(plyr)
library(MCMCglmm)

list.files(path="./data/delph")
Gdat = read.table("./data/delph/traitdata.txt", header=T)
crossID = read.table("./data/delph/crossID.txt", header=T)

# Format the trait data
calyx_w3_F = rep(NA, nrow(Gdat))
calyx_w3_M = rep(NA, nrow(Gdat))
calyx_w4_F = rep(NA, nrow(Gdat))
calyx_w4_M = rep(NA, nrow(Gdat))
calyx_w5_F = rep(NA, nrow(Gdat))
calyx_w5_M = rep(NA, nrow(Gdat))

for(i in 1:nrow(Gdat)){
  if(Gdat$sex[i]=="F"){
    calyx_w3_F[i]=Gdat$calyx_w3[i]
    calyx_w4_F[i]=Gdat$calyx_w4[i]
    calyx_w5_F[i]=Gdat$calyx_w5[i]
  }
  else{
    calyx_w3_M[i]=Gdat$calyx_w3[i]
    calyx_w4_M[i]=Gdat$calyx_w4[i]
    calyx_w5_M[i]=Gdat$calyx_w5[i]
  }}

trdata = data.frame(calyx_w3_F,calyx_w3_M,calyx_w4_F,calyx_w4_M,calyx_w5_F,calyx_w5_M)
trdata$mean_calyx_F = rowMeans(subset(trdata, select=c("calyx_w3_F", "calyx_w4_F","calyx_w5_F")), na.rm=T)
trdata$mean_calyx_M = rowMeans(subset(trdata, select=c("calyx_w3_M", "calyx_w4_M","calyx_w5_M")), na.rm=T)

data = data.frame(Gdat[,1:3], trdata)
names(data)[1] = "Cross"
data = data[-1020,] #Remove duplicated ID

crossID = crossID[,1:3]
data = merge(data, crossID, by="Cross", all.x=T)
data$animal = paste(data$Cross, data$individual, sep="-")

head(data)

# Format the pedigree data
pedigree = subset(data, select=c("animal","Female","Male"))
names(pedigree) = c("animal", "dam", "sire")
parentped = data.frame(animal=unique(c(as.character(pedigree$dam), as.character(pedigree$sire))), dam=NA,sire=NA)
pedigree = rbind(parentped, pedigree)
pedigree$animal = factor(pedigree$animal)
pedigree$dam = factor(pedigree$dam)
pedigree$sire = factor(pedigree$sire)
head(pedigree, 30)

# Mean-scale and multiply by 10
head(data)
data[,c(4:11)]=apply(data[,c(4:11)], 2, function(x) 10*x/mean(x, na.rm=T))

data$dampop = factor(substr(data$Female, 1, 1))
data$sirepop = factor(substr(data$Male, 1, 1))
str(data)

# Run the MCMCglmm analysis
invA = inverseA(pedigree)$Ainv

# Two traits
n = 2
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*400
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G3=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

samples = 1000
thin = 50
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

a=Sys.time()
mod<-MCMCglmm(c(mean_calyx_F, mean_calyx_M) ~ -1+trait,
              random = ~us(trait):animal + us(trait):dampop + us(trait):sirepop,
              rcov = ~us(trait):units,
              data = data, ginverse = list(animal = invA),
              family = rep("gaussian", n), prior = prior, 
              nitt = nitt, burnin = burnin, thin = thin)
Sys.time()-a

save(mod, file="./analyses/delph_Silene/Gmat75k_2traits.RData")

load(file="./analyses/delph_Silene/Gmat75k_2traits.RData")
n=2

# Check convergence
summary(mod$VCV)
plot(mod$VCV[,4])

# Compile G matrix
gmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(gmat) = rownames(gmat) = c("calyx_F", "calyx_M")
round(gmat, 2)

# Compile D matrix (among dam pops)
dmat = matrix(apply(mod$VCV, 2, median)[5:8], nrow=n)
colnames(dmat) = rownames(dmat) = c("calyx_F", "calyx_M")
dmat

# Compile S matrix (among sire pops)
smat = matrix(apply(mod$VCV, 2, median)[9:12], nrow=n)
colnames(smat) = rownames(smat) = c("calyx_F", "calyx_M")
smat

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)
evolvabilityMeans(smat)

round(gmat, 2)
round(dmat, 2)
round(smat, 2)

# Compute eigenvectors etc.
g_ev = eigen(gmat)$vectors
var_g_g = evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
var_s_g = evolvabilityBeta(smat, Beta = g_ev)$e

d_ev = eigen(dmat)$vectors
var_g_d = evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d = evolvabilityBeta(dmat, Beta = d_ev)$e
var_s_d = evolvabilityBeta(smat, Beta = d_ev)$e

s_ev = eigen(smat)$vectors
var_g_s = evolvabilityBeta(gmat, Beta = s_ev)$e
var_d_s = evolvabilityBeta(dmat, Beta = s_ev)$e
var_s_s = evolvabilityBeta(smat, Beta = s_ev)$e

p_ev = eigen(pmat)$vectors
var_g_p = evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p = evolvabilityBeta(dmat, Beta = p_ev)$e

# Compute summary stats
mg = lm(log(var_d_g)~log(var_g_g))
beta_g = summary(mg)$coef[2,1]
beta_g
#r2_g = summary(mg)$r.squared
#r2_g

md = lm(log(var_d_d)~log(var_g_d))
beta_d = summary(md)$coef[2,1]
beta_d
#r2_d = summary(md)$r.squared
#r2_d

ms = lm(log(var_s_s)~log(var_g_s))
beta_s = summary(ms)$coef[2,1]
beta_s
#r2_d = summary(md)$r.squared
#r2_d

mp = lm(log(var_d_p)~log(var_g_p))
beta_p = summary(mp)$coef[2,1]
beta_p
r2_p = summary(mp)$r.squared
r2_p

# Plot (dmat)
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Silene latifolia", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
#points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

# Angles
acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

# Plot (smat)
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_s), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_s), na.rm=T))
ymin = log10(min(c(var_s_g, var_s_s), na.rm=T))
ymax = log10(max(c(var_s_g, var_s_s), na.rm=T))
plot(log10(var_g_g), log10(var_s_g), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [x100])", 
     main="Silene latifolia", las=1)
points(log10(var_g_s), log10(var_s_s), pch=16)
#points(log10(var_g_p), log10(var_d_p), pch=16, col="green")
points(log10(diag(gmat)), log10(diag(smat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

# Angles
acos(t(g_ev[,1]) %*% s_ev[,1])*(180/pi)

#### Divergence vectors ####

means = apply(data[, 10:11], 2, function(x) tapply(x, data$dampop, mean, na.rm=T))
means
z0 = colMeans(means)

# Compute divergence etc.
outdat=matrix(NA, nrow=nrow(means), ncol=3)
for(i in 1:nrow(means)){
  z1=unlist(means[i,])
  delta = log(z1)-log(z0)
  scale_delta = delta/sqrt(sum(delta^2)) 
  d = delta
  div = mean(abs(d))*100
  e_delta = evolvabilityBeta(gmat, scale_delta)$e
  c_delta = evolvabilityBeta(gmat, scale_delta)$c
  
  outdat[i,]=c(div, e_delta, c_delta)
}

data.frame(means, outdat)

x11(width=5, height=5)
par(mar=c(4,4,5,4))
plot(outdat[,1], outdat[,2], pch=16, ylim=c(0,2), las=1,
     xlab="", ylab="",
     main="Silene latifolia")
mtext("Divergence from focal population (x100)", 1, line=2.5)
mtext("Evolvability (%)", 2, line=2.5)
points(outdat[,1], outdat[,3], pch=16, col="grey")

evolvabilityMeans(gmat)
abline(h=evolvabilityMeans(gmat)[1], lty=2)
abline(h=evolvabilityMeans(gmat)[4], lty=2, col="grey")
abline(h=evolvabilityMeans(gmat)[2], lty=1)
abline(h=evolvabilityMeans(gmat)[3], lty=1)
x = 5.5
text(x=x, y=evolvabilityMeans(gmat)[1], labels=expression(bar(e)), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[2], labels=expression(e[min]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[3], labels=expression(e[max]), xpd=T, cex=.8, adj=0)
text(x=x, y=evolvabilityMeans(gmat)[4], labels=expression(bar(c)), xpd=T, cex=.8, adj=0)
legend("topleft", c("e","c"), pch=16, col=c("black","grey"), bty="n")


#### Raw D matrix ####
head(data)
unique(data$Female)
data$dampop=substr(data$Female,1,1)
data$sirepop=substr(data$Male,1,1)

reddat=data[which(data$dampop==data$sirepop),]
head(reddat)
popmeans=ddply(reddat,.(dampop, sirepop), summarize,
               mean_calyx_F=mean(mean_calyx_F, na.rm=T),
               mean_calyx_M=mean(mean_calyx_M, na.rm=T))

popmeans

Dmat=cov(log(popmeans[,3:4]))*100
Dmat

