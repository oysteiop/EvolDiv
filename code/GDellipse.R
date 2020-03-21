GDellipse = function(dmat, gmat, location=c(0,0,0,0), ...){
  ctr    <- location[1:2]  # data centroid
  angles <- seq(0, 2*pi, length.out=500) 
  A = dmat
  eigVal  <- eigen(A)$values
  eigVec  <- eigen(A)$vectors
  eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
  #eigScl  <- eigVec %*% diag(eigVal)  # scale eigenvectors to length = square-root
  xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
  yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
  ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
  #ellBase <- cbind(eigVal[1]*cos(angles), eigVal[2]*sin(angles)) # normal ellipse
  ellRot  <- eigVec[,1:2] %*% t(ellBase)                                          # rotated ellipse
  
  x11(width=5, height=5)
  plot((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, 
       xlab="Trait 1", ylab="Trait 2",
       type="l", lwd=1, col="lightblue", las=1)
  polygon((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, col="lightblue")
  matlines(xMat[,1:2], yMat[,1:2], lty=1, lwd=1, col="black")
  
  A = gmat
  ctr = location[3:4]
  eigVal  <- eigen(A)$values
  eigVec  <- eigen(A)$vectors
  eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
  #eigScl  <- eigVec %*% diag(eigVal)  # scale eigenvectors to length = square-root
  xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
  yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
  ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
  #ellBase <- cbind(eigVal[1]*cos(angles), eigVal[2]*sin(angles)) # normal ellipse
  ellRot  <- eigVec[,1:2] %*% t(ellBase)                                          # rotated ellipse
  polygon((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, lwd=1, col=rgb(1,0,0,.5))
  matlines(xMat[,1:2], yMat[,1:2], lty=1, lwd=1, col="black")
}
