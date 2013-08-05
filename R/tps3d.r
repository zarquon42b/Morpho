tps3d<-function(M,refmat,tarmat,lambda=0)
{   
    q<-dim(M)[1]
    p<-dim(refmat)[1]
    m<-dim(refmat)[2]
    if (m==3)
    {
    Lall<-CreateL(refmat,lambda=lambda)
    }
    else
    {Lall<-CreateL2D(refmat)
    }
    Linv<-Lall$Linv
    m2<-rbind(tarmat,matrix(0,m+1,m))
    coeff<-matrix(NA,p+m+1,m)
    transM<-matrix(NA,q,m)
    coeff<-Linv%*%m2
    transM<-.fx(refmat,M,coeff)
 
    return(transM)
    
}

  
