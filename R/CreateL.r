CreateL<-function(matrix,lambda=0)
{  
	 k<-dim(matrix)[1]
    q1<-matrix(c(rep(1,k)),k,1)
    K<-matrix(0,k,k)
    Q<-cbind(q1,matrix)
    O<-matrix(c(rep(0,16)),4,4)

storage.mode(K)<-"double"
storage.mode(matrix)<-"double"
    K<-.Fortran("createL",K,nrow(K),matrix,ncol(matrix))[[1]]
	
	diag(K)<-lambda
    	L<-rbind(cbind(K,Q),cbind(t(Q),O))
	
	L1<-try(solve(L),silent=TRUE)
    	if (class(L1)=="try-error")
		{cat("CreateL: singular matrix: general inverse will be used.\n")
		L1<-mpinv(L)		
		}

	Lsubk<-L1[1:k,1:k]

    Lsubk3<-rbind(cbind(Lsubk,matrix(0,k,2*k)),cbind(matrix(0,k,k),Lsubk,matrix(0,k,k)),cbind(matrix(0,k,2*k),Lsubk))
    return(list(L=L,Linv=L1,Lsubk=Lsubk,Lsubk3=Lsubk3))
}
