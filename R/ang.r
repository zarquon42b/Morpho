ang <- function(x,y) {
    if (!is.matrix(x) || !is.numeric(x))
        stop("x must be a numeric matrix")
    if (!is.numeric(y) || length(y) != ncol(x))
        stop("y must be a vector of ncol(x)")
    a <- .Call("ang_calcC",x,y)
    return(a)
}
angM <- function(x,y) {
    if (!is.matrix(x) || !is.numeric(x))
        stop("x must be a numeric matrix")
    if (!is.matrix(y) || ncol(y) != ncol(x))
        stop("y must be a matrix of same dimensionality as x")
    a <- .Call("ang_calcM",x,y)
    return(a)
}
