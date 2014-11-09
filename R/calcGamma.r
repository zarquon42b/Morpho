#' @importFrom Matrix forceSymmetric tcrossprod
calcGamma <- function(Gamma0,Lsubk3,U,dims,stepsize=1)
  {
      
      U <- as(U,"sparseMatrix")
      tUL <- crossprod(U,Lsubk3)
      ULU <- forceSymmetric(tUL%*%U)
      B <- tUL%*%Gamma0
      B <- as(B,"sparseMatrix")
      T <- solve(ULU,B)
      ULUT <- U%*%T
      Gamma0 <- Gamma0-stepsize*ULUT
      Gamma0 <- matrix(Gamma0,length(Gamma0)/dims,dims)
      return(Gamma0)
  }

calcProcDGamma <- function(U,Gamma0,mshape,dims,stepsize=1) {
    Tpart <- tcrossprod(U)
    mshape <- as.vector(mshape)
    tmpdiff <- Gamma0-mshape
    slided <- Gamma0-stepsize*(Tpart%*%tmpdiff)
    slided <- matrix(slided,length(slided)/dims,dims)
    return(slided)
}
