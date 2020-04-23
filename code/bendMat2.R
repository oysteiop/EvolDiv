bendMat2 = function(x, value = 0.01){
  e = eigen(x)
  L = diag(e$values[which(e$values>0)])
  newmat = e$vectors[,1:ncol(L)] %*% L %*% t(e$vectors[,1:ncol(L)]) 
  return(newmat)
}
