###########################################################
#### - Attenuation bias for univariate meta-analysis - ####
###########################################################
rm(list=ls())

library(devtools)
library(plyr)
library(reshape2)
library(lme4)
library(glmmTMB)
library(MCMCglmm)
library(evolvability)
library(MuMIn)
library(fBasics)

# Read datafiles
ddat = read.table("data/dmatdata.txt", header=T)
ddat$ID = paste(ddat$reference, ddat$species, ddat$environment, sep="_")
maxdists = read.csv(file="data/maxdists.csv")

# List of studies
studies = sort(unique(ddat$ID))
splist = tapply(as.character(ddat$species), ddat$ID, function(x) x[1])

# Compute among-pop variance for each trait in each study
outlist = list()
for(s in 1:length(studies)){
  red = ddat[ddat$ID==studies[[s]],]
  red$trait = factor(red$trait)
  df=data.frame(study_ID = studies[[s]], 
                species = splist[[s]], 
                trait = sort(unique(red$trait)),
                environment = sort(unique(red$environment)),
                npop = tapply(red$mean>-Inf, red$trait, sum, na.rm=T),
                d = cbind(tapply(log(abs(red$mean)), red$trait, var, na.rm=T)),
                mean = cbind(tapply(red$mean, red$trait, mean, na.rm=T)),
                ref = unique(red$reference))
  df = na.omit(df)
  outlist[[s]] = df
}
ddf = rbind.fill(outlist)

ddat$species_trait = paste0(ddat$species, "_", ddat$trait)
ddf$species_trait = paste0(ddf$species, "_", ddf$trait)

# Combine with data from Evolvability database
edat = read.table("data/evolvabilitydatabase2020.txt", header=T)
edat$species_measurement = paste0(edat$species,"_", edat$measurement)

evals = NULL
n_e = NULL
evar = NULL
for(i in 1:nrow(ddf)){
  w = which(as.character(edat$species)==as.character(ddf$species)[i] & as.character(edat$species_measurement)==as.character(ddf$species_trait)[i])
  evals[i] = mean(edat$evolvability[w], na.rm=T)
  n_e[i] = sum(!is.na(edat$evolvability[w]), na.rm=T)
  evar[i] = var(log(edat$evolvability[w]), na.rm=T)
}
ddf$evals = evals

sum(n_e>1)
sum(n_e>0)
ddf$var.obs = evar # Variance in log evolvability among repeated estimates
ddf$log_e = log(evals) # Log evolvability
ddf$n_e = n_e # Number of evolvability estimates
pop.mean.obs <- weighted.mean(x = ddf$var.obs,
                              w = n_e/sum(n_e, na.rm = T), 
                              na.rm = T)

ddf = ddf[-which(is.infinite(ddf$log_e)),] #Removing one -Inf value

pop.rel.k = 1 - sum(((ddf$log_e^2)*ddf$var.obs/(var(ddf$log_e, na.rm=T)-
                  pop.mean.obs+ddf$var.obs)), na.rm = T)/sum(ddf$log_e^2, na.rm = T)
pop.rel.k
