.fx <- function(refmat,M,coefs,time=TRUE,threads=1) { 	
    M <- cbind(1,M)
    splM <- .Call("tpsfx",refmat,  M, t(coefs),threads)
    return(splM)
}
