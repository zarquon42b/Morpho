calcGamma <- function(Gamma0,Lsubk3,U,dims,weights=NULL)
    {
        if (require(Matrix)) {
            U <- as( U,"sparseMatrix")
            ULU <- Matrix::crossprod(U,(Lsubk3%*%U))
            ULU2 <-  Matrix::crossprod(U,(Lsubk3%*%Gamma0))
            diag(ULU) <- Matrix::diag(ULU)+1e-8
            if (is.null(weights))
                weights <- 1
            Gamma1 <- try(Gamma0 - weights*(U%*%Matrix::solve(ULU,ULU2)))
            if (class(Gamma1)=="try-error") {
                cat("calcGamma: singular matrix: general inverse will be used.\n")
                B <- armaGinv(as.matrix(ULU))		
                Gamma1 <- Gamma0-weights*(U%*%B%*%crossprod(U,(Lsubk3%*%Gamma0)))
            }
            Gamatrix <- matrix(Gamma1,length(Gamma1)/dims,dims)
            return(list(Gamma1=Gamma1,Gamatrix=Gamatrix))
        } else
            stop("please install package Matrix")
    }
