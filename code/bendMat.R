bendMat = function(x, value = 0.01){
  e = eigen(x)
  L = diag(e$values)
  diag(L)[which(diag(L)<0)] = value
  newmat = e$vectors %*% L %*% t(e$vectors) 
  return(newmat)
}
