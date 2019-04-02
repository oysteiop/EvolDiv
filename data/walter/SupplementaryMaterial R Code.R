
library(MASS);library(MCMCglmm);library(gdata);library(matrixcalc)

#### Equation 3 - for Dune ecotype only ####
dat <- read.csv("EcotypeData.csv", header = TRUE) 
names(dat)

dat <- droplevels(subset(dat, Type == "Dune"))

Dune <- data.frame(Block = dat$Block, Type = dat$Type, sire = dat$sire, dam = dat$dam, apply(dat[,8:17], 2, function(x) scale(x, center = TRUE, scale = TRUE)))

p <- diag(10)

p1 <- list(R = list(V = p*0.25, nu = 9.002),
           G = list(G1 = list(V = diag(10), nu = 10, alpha.mu = rep(0, 10), alpha.V = p),
                    G2 = list(V = diag(10), nu = 10, alpha.mu = rep(0, 10), alpha.V = p),
                    G3 = list(V = diag(10), nu = 10, alpha.mu = rep(0, 10), alpha.V = p)))

model_Dune <- MCMCglmm(cbind(Height, MSL_W, SB, MSD, Area, P2A2, Circularity, Nindents.Peri, IndentWidthMean, IndentDepthMean) ~ trait * Block,
                   random = ~ us(trait):sire + us(trait):dam + us(trait):sire:dam ,
                   rcov = ~ us(trait):units,
                   family = rep("gaussian",10),
                   pr = T,
                   nitt = 2100000, thin = 2000, burnin = 100000,
                   prior = p1,
                   data = Dune)
save(model_Dune, file = "model_Dune.Rdata")

#### Construct matrices ####
load("model_Dune.Rdata")
load("model_Head.Rdata")
load("model_Table.Rdata")
load("model_Wood.Rdata")

EcoS <- array(,c(10,10,4,1000))
for (i in 1:1000){
  EcoS[,,1,i] <- matrix(model_Dune$VCV[i,1:100], ncol = 10)
  EcoS[,,2,i] <- matrix(model_Head$VCV[i,1:100], ncol = 10)
  EcoS[,,3,i] <- matrix(model_Table$VCV[i,1:100], ncol = 10)
  EcoS[,,4,i] <- matrix(model_Wood$VCV[i,1:100], ncol = 10)
}

EcoG <- EcoS * 4

#### Covariance Tensor - For full tutorial see supplementary material in Aguirre et al. (2014b) ####
n = 10
m = 4
MCMCsamp = 1000
Gnames <- c("Dune","Headland","Tableland","Woodland")
traitnames <- c("Height", "MSL_W", "SB", "MSD",  "Area", "P2A2", "Circularity", "Nindents", "IndentWidthMean", "IndentDepthMean")

#START
covtensor <- function(Gs){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  neigten <- n*(n+1)/2 
  #Number of eigentensors
  MCMC.S <- array(,c(neigten, neigten, MCMCsamp))
  dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
  for (k in 1:MCMCsamp){
    MCMCG <- Gs[,,,k] 
    MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
    #find the variances of the kth G and store them 
    MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
    #find the covariances of the kth G and store them
    MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
    #fill the upper left quadrant of the kth S
    MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
    #fill the lower right quadrant of the kth S
    MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
    #fill the upper right quadrant of the kth S
    MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
    #fill the lower left quadrant of the kthS
  }  
  av.S <- apply(MCMC.S, 1:2, mean)
  #posterior mean S
  av.S.val <- eigen(av.S)$values
  #eigenvalues of posterior mean S 
  av.S.vec <- eigen(av.S)$vectors
  #eigenvalues of posterior mean S
  eTmat <- array(, c(n, n, neigten))
  dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
  for (i in 1:neigten){
    emat <- matrix(0, n, n) 
    lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- av.S.vec[1:n,i]
    eTmat[,,i] <- emat 
  }
  #construct the second-order eigentensors of posterior mean S
  eT.eigen <- array(, c(n+1, n, neigten))
  for (i in 1:neigten){
    eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
    #Eigenvalues of the ith eigentensor
    eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
    #Eigenvectors of the ith eigentensor
    eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = TRUE), i]
  }
  MCMC.S.val <- matrix(, MCMCsamp, neigten)
  colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
  for (i in 1:MCMCsamp){
    for(j in 1:neigten){
      MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
    }
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
  av.G.coord <- array(, c(m, neigten, 1))
  dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
  tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
  colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
  rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}
#END

nnonzero <- min(n*(n+1)/2,m-1)

Gtensor <- covtensor(EcoG)

Gtensor$tensor.summary #### Tensor summary
Gtensor$MCMC.S.val #### The S matrix
Gtensor$MCMC.G.coord #### The coordinates of each matrix in each eigentensor

