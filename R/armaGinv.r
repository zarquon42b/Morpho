#' calculate Pseudo-inverse of a Matrix using RcppArmadillo
#'
#' a simple wrapper to call Armadillo's pinv function
#' @param x numeric matrix
#' @param tol numeric: maximum singular value to be considered
#' @return Pseudo-inverse
#' @examples
#' mat <- matrix(rnorm(12),3,4)
#' pinvmat <- armaGinv(mat)
#' @export
armaGinv <- function(x, tol=NULL)
    {
        if (!is.matrix(x) || !is.numeric(x))
            stop("input must be a matrix")
        out <- .Call("armaGinv", x ,tol)
        return(out)
    }
