qqmat <- function(x,output=FALSE)
  {
    x <- as.matrix(x)
    center <- colMeans(x)
    n <- nrow(x); p <- ncol(x); cov <- cov(x);
    d <- mahalanobis(x,center,cov) # distances
    x1 <- qchisq(ppoints(n),df=p)
    xlab <- "expected values for multivariate normal distribution"
    out <- qqplot(x1,d,ylim = c(range(x1)),main="QQ Plot Assessing Multivariate Normality",ylab="Mahalanobis D2",xlab=xlab)
 abline(0,1)
    if(output)
      {
        return(out)
  }
}
