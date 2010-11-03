calcGamma<-function(Gamma0,Lsubk3,U,dims)
       {
	ULU<-crossprod(U,(Lsubk3%*%U))
	dU<-dim(ULU)[1]
	dia0<-as.matrix(diag(c(rep(1e-8,dU))))
	ULU<-ULU+dia0
	#B<-solve(ULU)
	B<-try(solve(ULU),silent=TRUE)
    	if (class(B)=="try-error")
		{cat("calcGamma: singular matrix: general inverse will be used.\n")
		B<-mpinv(ULU)		
		}
	Gamma1<-Gamma0-U%*%B%*%crossprod(U,(Lsubk3%*%Gamma0))
        Gamatrix<-matrix(Gamma1,length(Gamma1)/dims,dims)
        return(list(Gamma1=Gamma1,Gamatrix=Gamatrix))
        }
