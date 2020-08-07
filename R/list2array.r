#' converts a list of matrices to an array
#'
#' converts a list of matrices to an array
#' @param x a list containing matrices of the same dimensionality
#' @return returns an array concatenating all matrices
#' @export
list2array <- function(x) {
    xclass <- sapply(x,inherits,"matrix")
    classchk <- prod(xclass)
    if (!classchk)
        stop("all list entries must be matrices")
    xdim <- sapply(x,dim)
    dimchk <- prod(xdim[1,] == xdim[1,1])*prod(xdim[2,] == xdim[2,1])
    if (!dimchk)
        stop("all list entries must have the same dimensions")
    
    arr <- array(0,dim=c(dim(x[[1]]),length(x)))
    dimnames(arr)[[3]] <- names(x)
    if (!is.null(rownames(x[[1]])))
        dimnames(arr)[[1]] <- rownames(x[[1]])
     if (!is.null(colnames(x[[1]])))
         dimnames(arr)[[2]] <- colnames(x[[1]])
    for (i in 1:length(x))
        arr[,,i] <- x[[i]]
    return(arr)
}

#' reverts list2array, converting an array to a list of matrices 
#'
#' reverts list2array, converting an array to a list of matrices 
#' @param x array 
#' @return returns a list containing the matrices
#' @export
array2list <- function(x) {
    outlist <- list()
    for (i in 1:dim(x)[3]) {
        outlist[[i]] <- x[,,i]
    }
    names(outlist) <- dimnames(x)[[3]]
    return(outlist)
        
    
}
