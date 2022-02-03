bendMat = function(x, value=NULL){
  e = eigen(x)
  L = diag(e$values)
  if(is.null(value)){value = min(diag(L)[which(diag(L)>0)])/2}
  diag(L)[which(diag(L)<0)] = value
  newmat = e$vectors %*% L %*% t(e$vectors) 
  return(newmat)
}
