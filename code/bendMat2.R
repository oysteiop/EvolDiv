g = out$G
d = out$D

alignMat = function(g, d){
  ge <- eigen(g)
  G_ge <- diag(ge$values[which(ge$values>0)]) #Diagonal matrix with the positive G eigenvalues
  D_ge <- (t(ge$vectors)%*% d %*% ge$vectors)[1:ncol(G_ge), 1:ncol(G_ge)] #D matrix rotated along the G eigenvectors
  
  de <- eigen(D_ge)
  D_de <- diag(de$values[which(de$values>0)]) #Diagonal matrix with the positive D eigenvalues
  G_de <- (t(de$vectors) %*% G_ge %*% de$vectors)[1:ncol(D_de), 1:ncol(D_de)] #G matrix rotated along the D eigenvectors
  
  return(list(G_ge=G_ge, D_ge=D_ge, D_de=D_de, G_de=G_de))
}

a = alignMat(out$G, out$D)


alignMat2 = function(g, d){
  ge <- eigen(g)
  G_ge <- diag(ge$values[which(ge$values>0)]) #Diagonal matrix with the positive G eigenvalues
  D_ge <- (t(ge$vectors)%*% d %*% ge$vectors)[1:ncol(G_ge), 1:ncol(G_ge)] #D matrix rotated along the G eigenvectors
  
  de <- eigen(D_ge)
  D_de <- diag(de$values[which(de$values>0)]) #Diagonal matrix with the positive D eigenvalues
  G_de <- (t(de$vectors) %*% G_ge %*% de$vectors)[1:ncol(D_de), 1:ncol(D_de)] #G matrix rotated along the D eigenvectors
  
  ge <- eigen(G_de)
  G_ge <- diag(ge$values) 
  D_ge <- (t(ge$vectors)%*% D_de %*% ge$vectors) #D matrix rotated along the G eigenvectors
  
  return(list(G_ge=G_ge, D_ge=D_ge, D_de=D_de, G_de=G_de))
}

a2=alignMat2(out$G, out$D)

a$G_ge
a2$G_ge

a$D_ge
a2$D_ge

a$G_de
a2$G_de

a$D_de
a2$D_de

a
a2

plot(diag(a$G_de)[1:6], diag(a2$G_de))

a=a2
# Scaling relationship for the original traits
summary(lm(log(diag(out$D)) ~ log(diag(out$G))))$coef[2,1]

# Scaling relationship for the G eigenvectors
summary(lm(log(diag(a$D_ge)) ~ log(diag(a$G_ge))))$coef[2,1]

# Scaling relationship for the D eigenvectors
summary(lm(log(diag(a$D_de)) ~ log(diag(a$G_de))))$coef[2,1]

# D eigenvectors, conditional evolvability
G_de_c = evolvabilityBeta(a$G_de, eigen(a$D_de)$vectors)$c
summary(lm(log(diag(a$D_de)) ~ log(G_de_c)))$coef[2,1]






ge <- eigen(out$G)

reducedG <- diag(ge$values[which(ge$values>0)]) #Diagonal matrix with G eigenvalues

reducedD <- (t(ge$vectors)%*%out$D%*% ge$vectors)[1:ncol(reducedG), 1:ncol(reducedG)] #Added the t() and changed newG to reducedG


plot(log(diag(reducedG)), log(diag(reducedD)))
summary(lm(log(diag(reducedD))~log(diag(reducedG))))


#Her er diagonalane i reducedG og reducedD langs eigenvektorane til G, du kan rotere desse 
#dimensionsreduserte matrisene slik at dei er langs eigenvektorane til D ved

de <- eigen(reducedD)
D_de <- diag(de$values[which(de$values>0)])
G_de <- (t(de$vectors)%*% reducedG%*% de$vectors)[1:ncol(D_de), 1:ncol(D_de)]

summary(lm(log(diag(D_de))~log(diag(G_de))))


