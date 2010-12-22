mc.tps3d<-function(M,refmat,tarmat,lambda=0)
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
   # Linv<-try(solve(L),silent=TRUE)
#	if (is.character(Linv))
#		{cat("singular matrix: general inverse will be used.\n")
#		Linv<-mpinv(L)		
#		}
    m2<-rbind(tarmat,matrix(0,m+1,m))
    coeff<-matrix(NA,p+m+1,m)
    transM<-matrix(NA,q,m)
    for (i in 1:m)
      {
        coeff[,i]<-Linv%*%m2[,i]

      }
	#cat("calculating x displacement\n")
    	a.list<-as.list(1:m)
	warpM<-function(i)
		{
		transo<-mc.fx(refmat,M,coeff[,i])
		return(transo)
		}
	#if (m == 3)
	#	{dimo<-c("y","z")
	#	}
	#else
	#	{dimo<-c("y")
	#	}
	#for (i in 2:m)
    #{	cat(paste("calculating",dimo[i-1],"displacement\n"))
    #    transM[,i]<-mc.fx(refmat,M,coeff[,i],time=FALSE)
    #}
	a.list<-mclapply(a.list,warpM)
	for(i in 1:m)
	{transM[,i]<-a.list[[i]]}

    return(transM)
    
}
  
