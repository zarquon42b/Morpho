tps3d<-mc.tps3d<-function(M,refmat,tarmat)
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
    Linv<-Lall$Linv
   # Linv<-solve(L)
    m2<-rbind(tarmat,matrix(0,m+1,m))
    coeff<-matrix(NA,p+m+1,m)
    transM<-matrix(NA,q,m)
    #for (i in 1:m)
    #  {
        coeff<-Linv%*%m2

     # }
    #for (i in 1:m)
    #{
        transM<-fx(refmat,M,coeff)
   # }
    return(transM)
    
}

  
