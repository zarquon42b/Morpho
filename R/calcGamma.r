#' @importFrom Matrix forceSymmetric
calcGamma <- function(Gamma0,Lsubk3,U,dims,weights=NULL)
  {
      U <- as(U,"sparseMatrix")
      ULU <- forceSymmetric(crossprod(U,(Lsubk3%*%U)))
      B <- (t(U))%*%Lsubk3%*%Gamma0
      B <- as(B,"sparseMatrix")
      T <- solve(ULU,B)
      ULUT <- U%*%T
      Gamma1 <- Gamma0-ULUT
      Gamatrix <- matrix(Gamma1,length(Gamma0)/dims,dims)
      return(list(Gamma1=Gamma1,Gamatrix=Gamatrix))
  }
