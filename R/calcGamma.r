#' @importFrom Matrix forceSymmetric
calcGamma <- function(Gamma0,Lsubk3,U,dims,weights=NULL)
  {
      U <- as(U,"sparseMatrix")
      tUL <- crossprod(U,Lsubk3)
      ULU <- forceSymmetric(tUL%*%U)
      B <- tUL%*%Gamma0
      B <- as(B,"sparseMatrix")
      T <- solve(ULU,B)
      ULUT <- U%*%T
      Gamma0 <- Gamma0-ULUT
      Gamma0 <- matrix(Gamma0,length(Gamma0)/dims,dims)
      return(Gamma0)
  }
