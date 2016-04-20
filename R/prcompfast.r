#' fast Principal Component Analysis (PCA)
#'
#' fast Principal Component Analysis (PCA)
#' @param x a numeric or complex matrix (or data frame) which provides
#'          the data for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should be returned
#' @param center a logical value indicating whether the variables should be shifted to be zero centered. Alternately, a vector of length
#' @param scale. a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. The default is \code{FALSE} for consistency with S, but in general scaling is advisable.  Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.  The          value is passed to \code{scale}. equal the number of columns of \code{x} can be supplied.  The value is passed to \code{scale}.
#' @param tol a value indicating the magnitude below which components should be omitted. (Components are omitted if their standard deviations are less than or equal to \code{tol} times the standard deviation of the first component.)  With the default null setting, no components are omitted.  Other settings for tol could be \code{tol = 0} or \code{tol = sqrt(.Machine$double.eps)}, which would omit essentially constant components.
#' @param ... arguments passed to or from other methods.
#' @note this function returns the same results as \code{prcomp} (apart from sign differences) but uses smarter matrix decompositions making it faster for nrow(x) >> ncol(x) and nrow(x) << ncol(x).
#' @return
#' \code{prcomp} returns a list with class \code{prcomp} containing the followin components:
#'
#'\item{sdev}{the standard deviations of the principal components (i.e.,
#'          the square roots of the eigenvalues of the
#'          covariance/correlation matrix, though the calculation is
#'          actually done with the singular values of the data matrix).}
#'
#'\item{rotation:}{ the matrix of variable loadings (i.e., a matrix whose columns
#'          contain the eigenvectors).  The function \code{princomp} returns
#'          this in the element \code{loadings}.}
#'
#'\item{x:}{ if \code{retx} is true the value of the rotated data (the centred
#'          (and scaled if requested) data multiplied by the \code{rotation}
#'          matrix) is returned.  Hence, \code{cov(x)} is the diagonal matrix
#'          \code{diag(sdev^2)}.  For the formula method, \code{napredict()} is
#'          applied to handle the treatment of values omitted by the
#'          \code{na.action}.}
#'
#'\item{center, scale:}{ the centering and scaling used, or \code{FALSE}}.
#' pcafast <- prcompfast(iris[,1:4])
#' pcadefault <- prcompfast(iris[,1:4])
#' ## check if both results are idential (ignoring the sign)
#' all.equal(lapply(pcafast,abs),lapply(pcadefault,abs))
#' @export
prcompfast <- function(x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, ...) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if (any(sc == 0)) 
        stop("cannot rescale a constant/zero column to unit variance")

    if(n < p) {
        XX <- tcrossprod(x)
        svdXX <- svd(XX,nu=0)
    
        usepc <- nonzero(svdXX$d,tol)
        Dinv <- rep(0,length(svdXX$d))
        Dinv[usepc] <- 1/sqrt(svdXX$d[usepc])
        pcaBasis <- as.matrix(t(x)%*%(svdXX$v%*%Matrix::Diagonal(length(Dinv),Dinv)))
        sdev <- sqrt(svdXX$d/(n-1))
    } else {
        XX <- crossprod(x)
        svdXX <- svd(XX,nv=0)
        usepc <- nonzero(svdXX$d,tol)
        pcaBasis <- svdXX$u
        sdev <- sqrt(svdXX$d/(n-1))
    }
    pcaBasis <- pcaBasis[,usepc]
    sdev <- sdev[usepc]
    
    dimnames(pcaBasis) <- list(colnames(x), paste0("PC", seq_len(ncol(pcaBasis))))
    r <- list(sdev = sdev, rotation = pcaBasis, center = if (is.null(cen)) FALSE else cen,  scale = if (is.null(sc)) FALSE else sc)
    if (retx) 
        r$x <- x %*% pcaBasis
    class(r) <- "prcomp"
    r
    }


nonzero <- function(d,tol=NULL) {
    rank <- length(d)
    if (!is.null(tol)) {
        rank <- sum(d > (d[1L] * tol))
    }
    return(1:rank)
}
