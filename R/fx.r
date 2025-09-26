.fx <- function(refmat,M,coefs,time=TRUE,tpskernel=0,threads=1) { 	
    M <- cbind(1,M)
    splM <- .Call("tpsfx",refmat,  M, t(coefs),tpskernel,threads)
    return(splM)
}
