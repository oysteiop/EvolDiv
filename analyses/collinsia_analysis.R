###########################################
#### Collinsia, data from Åsa Lankinen ####
###########################################

rm(list=ls())
library(evolvability)
library(plyr)
library(reshape2)
library(MCMCglmm)
library(lme4)
library(kinship2)

list.files(path="./data/lankinen")
dat1 = read.table("./data/lankinen/collinsia_z1.txt", header=T)[,4:7]
dat2 = read.table("./data/lankinen/collinsia_z2.txt", header=T)[,4:7]
dat3 = read.table("./data/lankinen/collinsia_z3.txt", header=T)[,4:7]
dat4 = read.table("./data/lankinen/collinsia_z4.txt", header=T)[,4:7]
dat1$Id = factor(paste(dat1$Population, dat1$Id, sep="_"))
dat2$Id = factor(paste(dat2$Population, dat2$Id, sep="_"))
dat3$Id = factor(paste(dat3$Population, dat3$Id, sep="_"))
dat4$Id = factor(paste(dat4$Population, dat4$Id, sep="_"))
dat1$Population = factor(dat1$Population)
dat2$Population = factor(dat2$Population)
dat3$Population = factor(dat3$Population)
dat4$Population = factor(dat4$Population)
dat1$Family = factor(paste(dat1$Population, dat1$Family, sep="_"))
dat2$Family = factor(paste(dat2$Population, dat2$Family, sep="_"))
dat3$Family = factor(paste(dat3$Population, dat3$Family, sep="_"))
dat4$Family = factor(paste(dat4$Population, dat4$Family, sep="_"))

dat = merge(dat1, dat2, by="Id", all=T)
dat = merge(dat, dat3, by="Id", all=T)
dat = merge(dat, dat4, by="Id", all=T)

names(dat)

dat$fam = apply(dat, 1, function(x) na.omit(unique(x[c(3,6,9,12)])))
dat$pop = apply(dat, 1, function(x) na.omit(unique(x[c(2,5,8,11)])))
dat$fam=factor(dat$fam)
dat$pop=factor(dat$pop)

dat=subset(dat, select=c("pop", "fam", "Id", 
                         "Day_of_flower_start", 
                         "Stage_of_anther_stigma_contact",
                         "Stage_of_stigma_receptivity",
                         "Keel_length_mm"))

head(dat)
str(dat)
summary(dat)

rd = droplevels(na.omit(dat[,c(1:4)]))
hist(rd[,4])
str(rd)
mod = lmer(log(Day_of_flower_start)~1 + (1|pop/fam), data = rd)
summary(mod)
4*0.03095/(0.10167 + 0.03095 + 0.07799) #H2 = 0.58
100*(4*0.03095) #e = 12.38

rd = droplevels(na.omit(dat[,c(1:3, 7)]))
hist(rd[,4])
str(rd)
mod = lmer(log(Keel_length_mm)~1 + (1|pop/fam), data = rd)
summary(mod)
4*0.002674/(0.002079 + 0.002674 + 0.007653) #H2 = 0.86
100*(4*0.002674) #e = 1.07

rd = droplevels(na.omit(dat[,c(1:3, 5)]))
hist(rd[,4])
str(rd)
mod = lmer(log(Stage_of_anther_stigma_contact)~1 + (1|pop/fam), data = rd)
summary(mod)
4*0.005321/(0.007931 + 0.005321 + 0.025961) #H2 = 0.54
100*(4*0.005321) #e = 2.13

plot(log10(c(12.38, 1.07, 2.13)), log10(c(10.167, 0.2079, 0.5321)))

# Summary stats
popmeans = ddply(dat, .(pop), summarize,
      DFSm = mean(Day_of_flower_start, na.rm=T),
      SASCm = mean(Stage_of_anther_stigma_contact, na.rm=T),
      SOSRm = mean(Stage_of_stigma_receptivity, na.rm=T),
      KLm = mean(Keel_length_mm, na.rm=T))
popmeans$pop = factor(popmeans$pop)
rownames(popmeans) = popmeans$pop
popmeans=popmeans[,-1]
write.csv2(melt(t(popmeans)), "popmeans.csv")

popvars = ddply(dat, .(pop), summarize,
                 DFSm = var(Day_of_flower_start, na.rm=T),
                 SASCm = var(Stage_of_anther_stigma_contact, na.rm=T),
                 SOSRm = var(Stage_of_stigma_receptivity, na.rm=T),
                 KLm = var(Keel_length_mm, na.rm=T))
popvars$pop = factor(popvars$pop)
rownames(popvars) = popvars$pop
popvars=popvars[,-1]
write.csv2(melt(t(popvars)), "popvars.csv")

popns = ddply(dat, .(pop), summarize,
                DFSm = sum(Day_of_flower_start>-Inf, na.rm=T),
                SASCm = sum(Stage_of_anther_stigma_contact>-Inf, na.rm=T),
                SOSRm = sum(Stage_of_stigma_receptivity>-Inf, na.rm=T),
                KLm = sum(Keel_length_mm>-Inf, na.rm=T))
popns$pop = factor(popns$pop)
rownames(popns) = popns$pop
popns=popns[,-1]
write.csv2(melt(t(popns)), "popns.csv")

reddat = dat[,-which(colnames(dat)=="Stage_of_stigma_receptivity")]
reddat = droplevels(na.omit(reddat))
str(reddat)

########################################################
###### Joint estimation of D and average G ######
########################################################
names(reddat)
Gdat = reddat
traits = c("Day_of_flower_start", "Stage_of_anther_stigma_contact", "Keel_length_mm")   

X = subset(Gdat, select=c(traits, "pop", "fam"))
X[,1:length(traits)] = apply(X[,1:length(traits)], 2, function(x) 10*log(x))
X = na.omit(X)
head(X)
str(X)

samples = 1000
thin = 100
burnin = samples*thin*.5
nitt = (samples*thin)+burnin

n = length(traits)
alpha.mu <- rep(0, n)
alpha.V <- diag(n)*50
prior<-list(R=list(V=diag(n), nu=n+0.002-1), 
            G=list(G1=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V),
                   G2=list(V=diag(n), nu=n, alpha.mu = alpha.mu, alpha.V = alpha.V)))

a = Sys.time()
mod<-MCMCglmm(cbind(Day_of_flower_start, Stage_of_anther_stigma_contact, Keel_length_mm) ~ -1+trait,
              random= ~us(trait):pop + us(trait):fam,
              rcov=~us(trait):units,
              data=X,
              family=rep("gaussian", n), prior=prior, nitt = nitt, burnin = burnin, thin = thin)
Sys.time() - a

save(mod, file="analyses/lankinen/mod150.RData")

plot(mod$VCV[,1])

#The D matrix
Gdat = reddat
names(Gdat)

traits = c("Day_of_flower_start", "Stage_of_anther_stigma_contact", "Keel_length_mm")
n = length(traits)

load(file="analyses/lankinen/mod150.RData")

dmat = matrix(apply(mod$VCV, 2, median)[1:(n*n)], nrow=n)
colnames(dmat) = rownames(dmat) = traits
dmat

# The (within-pop mean) G matrix
gmat = matrix(apply(mod$VCV, 2, median)[((n*n)+1):((n*n)*2)], nrow=n)
colnames(gmat) = rownames(gmat) = traits
round(gmat, 2)
cov2cor(gmat)

# The (within-individual) residual matrix
rmat = matrix(apply(mod$VCV, 2, median)[(((n*n)*2)+1):((n*n)*3)], nrow=n)
colnames(rmat) = rownames(rmat) = traits
round(rmat, 2)

4*diag(gmat)/(diag(dmat) + diag(gmat))

pmatdat = na.omit(subset(Gdat, select=c(traits)))
pmat = cov(pmatdat)

evolvabilityMeans(gmat)
evolvabilityMeans(dmat)

#Compute eigenvectors etc.
g_ev=eigen(gmat)$vectors
var_g_g=evolvabilityBeta(gmat, Beta = g_ev)$e
var_d_g = evolvabilityBeta(dmat, Beta = g_ev)$e
#var_d_g = diag(t(g_ev) %*% dmat %*% g_ev)

d_ev=eigen(dmat)$vectors
var_g_d=evolvabilityBeta(gmat, Beta = d_ev)$e
var_d_d=evolvabilityBeta(dmat, Beta = d_ev)$e
#var_d_d = diag(t(d_ev) %*% dmat %*% d_ev)

p_ev=eigen(pmat)$vectors
var_g_p=evolvabilityBeta(gmat, Beta = p_ev)$e
var_d_p=evolvabilityBeta(dmat, Beta = p_ev)$e

#Compute summary stats
mg=lm(log(var_d_g)~log(var_g_g))
beta_g=summary(mg)$coef[2,1]
beta_g
r2_g=summary(mg)$r.squared
r2_g

md=lm(log(var_d_d)~log(var_g_d))
beta_d=summary(md)$coef[2,1]
beta_d
r2_d=summary(md)$r.squared
r2_d

mp=lm(log(var_d_p)~log(var_g_p))
beta_p=summary(mp)$coef[2,1]
beta_p
r2_p=summary(mp)$r.squared
r2_p

#Plot both relationship,s upper/lower bounds of scaling relationship
#Plot
x11(width=5, height=5)
xmin = log10(min(c(var_g_g, var_g_d), na.rm=T))
xmax = log10(max(c(var_g_g, var_g_d), na.rm=T))
ymin = log10(min(c(var_d_g, var_d_d), na.rm=T))
ymax = log10(max(c(var_d_g, var_d_d), na.rm=T))
plot(log10(var_g_g), log10(var_d_g), 
     xlim=c(xmin, xmax), ylim=c(ymin-.5, ymax), 
     xlab="log10 (Evolvability [%])", 
     ylab="log10 (Divergence [%])", 
     main="Collinsia heterophylla", las=1)
points(log10(var_g_d), log10(var_d_d), pch=16)
points(log10(diag(gmat)), log10(diag(dmat)), pch=16, col="blue")
legend("bottomright", c("G eigenvectors", "D eigenvectors", "Traits"), pch=c(1,16, 16), col=c("black", "black", "blue"))

points(log10(var_g_p), log10(var_d_p), pch=16, col="red")

#Angles
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)
180-acos(t(g_ev[,1]) %*% d_ev[,1])*(180/pi)

# Ellipse plot
source("GDellipse.R")
GDellipse(dmat, gmat, xlim=c(-7,7), ylim=c(-7, 7), main="Senecio: Dune")
