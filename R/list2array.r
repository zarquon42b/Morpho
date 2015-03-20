#' converts a list of matrices to an array
#'
#' converts a list of matrices to an array
#' @param x a list containing matrices of the same dimensionality
#' @return returns an array concatenating all matrices
#' @export
list2array <- function(x) {
    xclass <- sapply(x,class)
    classchk <- prod(xclass == "matrix")
    if (!classchk)
        stop("all list entries must be matrices")
    xdim <- sapply(x,dim)
    dimchk <- prod(xdim[1,] == xdim[1,1])*prod(xdim[2,] == xdim[2,1])
    if (!dimchk)
        stop("all list entries must have the same dimensions")
    
    arr <- array(0,dim=c(dim(x[[1]]),length(x)))
    dimnames(arr)[[3]] <- names(x)
    for (i in 1:length(x))
        arr[,,i] <- x[[i]]
    return(arr)
}
