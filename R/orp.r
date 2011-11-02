orp<-function(A)
{ 
  p<-dim(A)[1]
  k<-dim(A)[2]
  n<-dim(A)[3]
  
  mshape1<-apply(A,c(1,2),mean)
  m.size <-c.size(mshape1)
  Y1<-as.vector(mshape1/m.size)
  oo<-as.matrix(rep(1,n))%*%as.vector(Y1)
  I<-diag(1,k*p)
  mat<-matrix(NA, n, k*p)
  for (i in 1:n)	{mat[i,]<-as.vector(A[,,i])/m.size}
  Xp<-mat%*%(I-(Y1%*%t(Y1)))
  Xp1<-Xp+oo
  return(proj=array(t(Xp1), dim=c(p, k, n)))
}
