orp <- function(A, mshape=NULL)
{ 
  p <- dim(A)[1]
  k <- dim(A)[2]
  n <- dim(A)[3]
  if (is.null(mshape))
      mshape <- arrMean3(A)
      
  m.size <- cSize(mshape)
  Xc <- as.vector(mshape/m.size)
  X <- vecx(A)
  ##direction along mshape onto plane
  XtoPlane <- t(apply(X,1,function(x){x <- t(crossprod(x,Xc)*Xc)}))
  X1 <- X-XtoPlane
  X1 <- t(X1)+Xc
  return(proj=array(X1, dim=c(p, k, n)))
}
