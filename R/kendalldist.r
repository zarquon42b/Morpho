kendalldist <- function(x,y)
    {
        m <- ncol(x)
        x <- apply(x,2,scale,scale=F)
        y <- apply(y,2,scale,scale=F)
        x <- x/cSize(x)
        y <- y/cSize(y)
        if (max(abs(x - y) > 0))
            {
                if (m == 3)
                    {
                        eigxy <- eigen(t(x)%*%tcrossprod(y)%*%(x))$values
                        signchk <- det(crossprod(y,x))
                        good <- which(eigxy > 0)
                        eigxy <- sqrt(eigxy[1:m])
                        eigxy[m] <- sign(signchk)*eigxy[m]
                        rho <- acos(min(sum(eigxy[1:m]),1))
                    }
                else
                    {
                        ## this is copied from 'riemdist' in shapes package
                        x <- realtocomplex(x)
                        y <- realtocomplex(y)
                        rho <- acos(min(1, (Mod(t(Conj(x)) %*% y))))
                    }
            }
                    else
            rho <- 0
            return(rho)
    }
