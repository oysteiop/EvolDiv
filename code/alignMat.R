alignMat = function(g, d){
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