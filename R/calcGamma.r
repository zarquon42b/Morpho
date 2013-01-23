calcGamma<-function(Gamma0,Lsubk3,U,dims,weights=NULL)
  {
    require(Matrix)
    U <- as( U,"sparseMatrix")
    ULU <- Matrix::crossprod(U,(Lsubk3%*%U))
                                        ##ULU <- t(U)%*%(Lsubk3%*%U)
	dU<-dim(ULU)[1]
	#dia0 <- (diag(c(rep(1e-8,dU))))
        #dia0 <- as(dia0,"sparseMatrix")
        dia0 <- sparseMatrix(1:dU,1:dU,x=1e-8,dims=c(dU,dU))
	ULU<-ULU+dia0
	#B<-solve(ULU)
	B<-try(solve(ULU),silent=TRUE)
    	if (class(B)=="try-error")
		{cat("calcGamma: singular matrix: general inverse will be used.\n")
		B<-ginv(ULU)		
               }
        B <- as(B,"sparseMatrix")
        if (is.null(weights))
          {
            weights <- 1
          }
        #tmp <- crossprod(U,(Lsubk3%*%Gamma0))
        #tmp <- crossprod(t(B),tmp)
        #tmp <- crossprod(t(U),tmp)
        #Gamma1 <- Gamma0-weights*tmp
        Gamma1<-Gamma0-weights*(U%*%B%*%Matrix::crossprod(U,(Lsubk3%*%Gamma0)))
        Gamatrix<-matrix(Gamma1,length(Gamma1)/dims,dims)
        return(list(Gamma1=Gamma1,Gamatrix=Gamatrix))
        }
