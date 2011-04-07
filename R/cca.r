cca <- function(X,Y,cor=T)
  {

    dim.x <- dim(X)
    dim.y <- dim(Y)
    part <- cbind(X,Y)
    if (cor)
      {
        
        co <- cor(part)
      }
    else
      {
        co <-  cov(part)
      }
    sxx <- co[1:dim.x[2],1:dim.x[2]]
    syy <- co[-c(1:dim.x[2]),-c(1:dim.x[2])]
    A <- solve(co[1:dim.x[2],1:dim.x[2]]) %*% co[1:dim.x[2],-c(1:dim.x[2])]
    B <- solve(co[-c(1:dim.x[2]),-c(1:dim.x[2])]) %*% co[-c(1:dim.x[2]),c(1:dim.x[2])]
    print(dim.x)
    corr <- Re(eigen(A%*%B)$values)
    corr <- sqrt(corr[which(corr > 1e-10)])
    print(corr)
    lambda <- det(co)/(det(sxx)*det(syy))
         print(lambda)              
    return(lambda)
  }
    
