.fx <- function(refmat,M,coefs,time=TRUE)
{ 	
    q <- dim(M)[1]
    p <- dim(refmat)[1]
    m <- dim(refmat)[2]
    M1 <- cbind(1,M)
    coefs <- t(coefs)
    storage.mode(M) <- "double"
    storage.mode(refmat) <- "double"
    storage.mode(coefs) <- "double"
                                        #splM <- .Fortran("tpsfx",refmat,p,M,q,M1,refmat[,1],coefs,M)[[8]]
    splM <- .Call("tpsfx",refmat, M, M1, coefs)
    
    return(splM)
}
