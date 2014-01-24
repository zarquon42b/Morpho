armaGinv <- function(x, tol= sqrt(.Machine$double.eps))
    {
        if (!is.matrix(x))
            stop("input must be a matrix")
        out <- .Call("armaGinv", x ,tol)
        return(out)
    }
