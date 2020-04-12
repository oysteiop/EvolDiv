Gvar2G = function(Gvar, Vp){
  vars = diag(Gvar)*Vp
  G = Gvar*sqrt(tcrossprod(vars))
  diag(G) = vars
  return(G)
}