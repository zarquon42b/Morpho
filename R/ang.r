ang <- function(x,y) {
    if (!is.matrix(x) || !is.numeric(x))
        stop("x must be a numeric matrix")
    if (!is.numeric(y) || length(y) != ncol(x))
        stop("y must be a vector of ncol(x)")
    a <- .Call("ang_calcC",x,y)
    return(a)
}
