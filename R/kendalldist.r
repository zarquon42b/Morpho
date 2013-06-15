kendalldist <- function(x,y)
    {
        m <- ncol(x)
        x <- apply(x,2,scale,scale=FALSE)
        y <- apply(y,2,scale,scale=FALSE)
        x <- x/cSize(x)
        y <- y/cSize(y)
        if (max(abs(x - y) > 0))
            {
                if (m == 3)
                    {
                        eigxy <- eigen(t(x)%*%tcrossprod(y)%*%(x),symmetric = TRUE)$values
                        signchk <- det(crossprod(y,x))
                        good <- which(eigxy > 0)
                        eigxy <- sqrt(eigxy[1:m])
                        eigxy[m] <- sign(signchk)*eigxy[m]
                        rho <- acos(min(sum(eigxy[1:m]),1))
                    }
                else
                    {
                        ## this is copied from 'riemdist' in shapes package
                        x <- x[,1]+(0+1i)*x[,2]
                        y <- y[,1]+(0+1i)*y[,2]
                        rho <- acos(min(1, (Mod(t(Conj(x)) %*% y))))
                    }
            }
                    else
            rho <- 0
            return(rho)
    }
