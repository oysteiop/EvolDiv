bendMat = function(x, value = 0.01){
  e = eigen(x)
  L = diag(e$values)
  value = min(diag(L)[which(diag(L)>0)])/10
  diag(L)[which(diag(L)<0)] = value
  newmat = e$vectors %*% L %*% t(e$vectors) 
  return(newmat)
}
