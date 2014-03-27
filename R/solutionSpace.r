#' returns the solution space (basis and translation vector) for an equation system
#'
#' returns the solution space (basis and translation vector) for an equation system
#'
#' @param A numeric matrix
#' @param b numeric vector
#'
#' @return
#' \item{basis}{matrix containing the basis of the solution space}
#' \item{translate}{translation vector}
#'
#' @details For a linear equationsystem, \eqn{Ax = b}{Ax = b}, the solution space then is
#' \deqn{x = A^* b + (I - A^* A) y}{x = A'b + (I - A' A)}
#' where \eqn{A^*}{A'} is the Moore-Penrose pseudoinverse of \eqn{A}{A}.
#' The QR decomposition of \eqn{I - A^* A}{I - A'A} determines the dimension of and basis of the solution space.
#'@examples
#' A <- matrix(rnorm(21),3,7)
#' b <- c(1,2,3)
#' subspace <- solutionSpace(A,b)
#' dims <- ncol(subspace$basis) # we now have a 4D solution space
#' ## now pick any vector from this space. E.g
#' y <- 1:dims
#' solution <- subspace$basis%*%y+subspace$translate # this is one solution for the equation above
#' A%*%solution ## pretty close
#' @export
solutionSpace <- function(A,b){
    Apinv <- armaGinv(A)
    translate <- c(Apinv%*%b)
    basis <- Apinv%*%A
    if (isTRUE(all.equal(diag(ncol(A)),basis))) {
        basis <- matrix(0,ncol(A),0)
    } else {
        basis <- -basis
        diag(basis) <- 1+diag(basis)
        qrcheck <- qr(basis)
        basis <- qr.Q(qrcheck)[,1:qrcheck$rank,drop=FALSE]
    }
    out <- list(basis=basis,translate=translate)
    return(out)
}
