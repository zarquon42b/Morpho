orp <- function(A, mshape=NULL)
{ 
  p <- dim(A)[1]
  k <- dim(A)[2]
  n <- dim(A)[3]
  if (is.null(mshape))
      mshape <- apply(A,c(1,2),mean)
      
  m.size  <- cSize(mshape)
  Xc <- as.vector(mshape/m.size)
  Ikp <- diag(k*p)
  X <- vecx(A)
  #X <- X/m.size ##remove resizing -> bad for large shape variation
  X1 <- X%*%(Ikp - tcrossprod(Xc))
  XcM <- as.matrix(rep(1,n))%*%Xc
  X1 <- X1+XcM
  return(proj=array(t(X1), dim=c(p, k, n)))
}
