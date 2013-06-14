scaleproc<-function (a3, proc.output = FALSE) 
{
#### this is a copy of the function "bgpa" from the shapes package ####
#### Copyright by Ian Dryden
  
    h <- 0
    zd <- a3
    s <- 0
    n <- dim(a3)[3]
    #traceXtX <- function(x) x <- sqrt(sum(diag(crossprod(x))))
    traceXtX <- function(x) x <- sum(x^2)
    aa <- apply(zd, c(3), traceXtX)
    s <- sum(aa)
    omat <- vecx(zd)
    kk <- dim(omat)[2]
    nn <- dim(omat)[1]
    if (nn > kk) {
        qq <- rep(0, times = nn)
        for (i in 1:n) {
            qq[i] <- var(omat[i, ]) * (n - 1)/n
            omat[i, ] <- omat[i, ] - mean(omat[i, ])
        }
        omat <- diag(sqrt(1/qq)) %*% omat
        n <- kk
        Lmat <- t(omat) %*% omat/n

	eig <- eigen(Lmat, symmetric = TRUE)
      
	 U <- eig$vectors
        lambda <- eig$values
        V <- omat %*% U
        vv <- rep(0, times = n)
        for (i in 1:n) {
            vv[i] <- sqrt(t(V[, i]) %*% V[, i])
            V[, i] <- V[, i]/vv[i]
        }
        delta <- sqrt(abs(lambda/n)) * vv
        od <- order(delta, decreasing = TRUE)
        delta <- delta[od]
        V <- V[, od]
        h <- sqrt(s/aa) * V[, 1]
    }
    if (kk >= nn) {
        zz <- cor(t(vecx(zd)))
        h <- sqrt(s/aa) * eigen(zz)$vectors[, 1]
    }
    h <- abs(h)

    return(h)
}

