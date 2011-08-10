pls_2B <- function(y,x)
  {
    if (length(dim(x)) == 3)
      {
        x <- vecx(x)
      }
    
    xdim <- dim(x)
    ydim <- dim(y)
    
    cova <- cov(cbind(x,y))
    svd.cova <- svd(cova[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])
    
    #fitted.values <- t(svd.cova$v%*%t(svd.cova$u)%*%t(x))
    z1 <- x%*%svd.cova$u
    fitted.values <- (fitted.values) %*%t(svd.cova$v)
    z2 <-  y%*%svd.cova$v
    cors <- 0
    compcov <- sum
    for(i in 1:length(svd.cova$d))
        {cors[i] <- cor(z1[,i],z2[,i])
       }
    
    out <- list(svd=svd.cova,z2=z2,z1=z1,cor=cors)
    return(out)
  }
