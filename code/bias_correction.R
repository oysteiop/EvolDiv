###########################################################
#### - Attenuation bias for univariate meta-analysis - ####
###########################################################
rm(list=ls())

library(plyr)
#library(reshape2)

# Read datafiles
ddat = read.csv2("data/dmatdata2.csv", dec=".", header=T)
ddat$ID = paste(ddat$reference, ddat$species, ddat$environment, sep="_")

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
edat = read.csv2("data/evolvabilitydatabase2020.csv", header=T, dec=".")
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

ddf$var.obs = evar # Variance in log evolvability among repeated estimates
ddf$log_e = log(evals) # Log evolvability
ddf$n_e = n_e # Number of evolvability estimates

ddf = ddf[-which(is.infinite(ddf$log_e)),] #Removing one -Inf value

ww = which(ddf$n_e>1)
ddf = ddf[ww,]
ddf = na.omit(ddf)

sum(ddf$n_e>1)

# Weighed observation error
pop.mean.obs <- weighted.mean(x = ddf$var.obs,
                              w = ddf$n_e/sum(ddf$n_e, na.rm = T), 
                              na.rm = T)

var.pred = 2.948813 # Variance of predictor in meta-analytical model


pop.rel.k = 1 - sum(((ddf$log_e^2)*ddf$var.obs/(var.pred-pop.mean.obs+ddf$var.obs)))/
  sum(ddf$log_e^2)

pop.rel.k
0.76/pop.rel.k # Corrected slope
0.76/0.777

# Simpler estimates
1-pop.mean.obs/var.pred
0.76/(1-pop.mean.obs/var.pred) # Corrected slope

1-mean(ddf$var.obs/var.pred)
0.76/(1-mean(ddf$var.obs/var.pred))
