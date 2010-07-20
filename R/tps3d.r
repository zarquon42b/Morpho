tps3d<-function(M,refmat,tarmat)
{   
    q<-dim(M)[1]
    p<-dim(refmat)[1]
    m<-dim(refmat)[2]
    if (m==3)
    {
    Lall<-CreateL(refmat)
    }
    else
    {Lall<-CreateL2D(refmat)
    }
    L<-Lall$L
    Linv<-solve(L)
    m2<-rbind(tarmat,matrix(0,m+1,m))
    coeff<-matrix(NA,p+m+1,m)
    transM<-matrix(NA,q,m)
    for (i in 1:m)
      {
        coeff[,i]<-Linv%*%m2[,i]

      }
    for (i in 1:m)
    {
        transM[,i]<-fx(refmat,M,coeff[,i])
    }
    return(transM)
    
}
  