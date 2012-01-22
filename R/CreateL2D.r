CreateL2D<-function(matrix)
{   k<-dim(matrix)[1]
    q1<-matrix(c(rep(1,k)),k,1)
    K<-matrix(0,k,k)
    Q<-cbind(1,matrix)
    O<-matrix(0,3,3)

    for (i in 1:k)
    {
      for (j in 1:k)
          {r2<-sum((matrix[i,]-matrix[j,])^2)
            K[i,j]<-r2*log(r2)   }
    }
    K[which(is.na(K))]<-0
    L<-rbind(cbind(K,Q),cbind(t(Q),O))
    
	L1<-try(solve(L),silent=TRUE)
    	if (class(L1)=="try-error")
		{cat("singular matrix: general inverse will be used.\n")
		L1<-mpinv(L)		
		}
    Lsubk<-L1[1:k,1:k]
    Lsubk3<-rbind(cbind(Lsubk,matrix(0,k,k)),cbind(matrix(0,k,k),Lsubk))
    return(list(L=L,Linv=L1,Lsubk=Lsubk,Lsubk3=Lsubk3))
}
