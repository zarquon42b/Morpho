armaGinv <- function(x, tol=NULL)
    {
        if (!is.matrix(x) || !is.numeric(x))
            stop("input must be a matrix")
        out <- .Call("armaGinv", x ,tol)
        return(out)
    }
