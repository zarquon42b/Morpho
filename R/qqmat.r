qqmat <- function(x,output=FALSE,square=TRUE)
  {
    x <- as.matrix(x)
    center <- colMeans(x)
    n <- nrow(x); p <- ncol(x); cov <- cov(x);
    d <- mahalanobis(x,center,cov) # distances
    x1 <- qchisq(ppoints(n),df=p)
  
    xlab <- "expected values for multivariate normal distribution"
    if (square)
      {
        ylim =  c(range(x1))
        out <- qqplot(x1,d,ylim =ylim ,main="QQ Plot Assessing Multivariate Normality",ylab="Mahalanobis D2",xlab=xlab)
      }
    else
      {
        out <- qqplot(x1,d ,main="QQ Plot Assessing Multivariate Normality",ylab="Mahalanobis D2",xlab=xlab)
      }
    out$d <- d
    abline(0,1)
    if(output)
      {
        return(out)
  }
}
