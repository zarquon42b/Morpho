#' Q-Q plot to assess normality of data
#' 
#' qqmat plots Mahalanobisdistances of a given sample against those expected
#' from a Gaussian distribution
#' 
#' 
#' @param x sample data: matrix or vector
#' @param output logical: if TRUE results are returned
#' @param square plot in a square window - outliers might be cut off.
#' @return if \code{output=TRUE}, the following values are returned
#' \item{x }{distances from an expected Gaussian distribution}
#' \item{y }{observed distances - sorted}
#' \item{d }{observed distances - unsorted}
#' @author Stefan Schlager
#' @seealso \code{\link{qqplot}}
#' 
#' @examples
#' 
#' require(MASS)
#' ### create normally distributed data
#' data <- mvrnorm(100,mu=rep(0,5),Sigma = diag(5:1))
#' qqmat(data)
#' 
#' ###create non normally distributed data
#' data1 <- rchisq(100,df=3)
#' qqmat(data1,square=FALSE)
#' 
#' @export
qqmat <- function(x,output=FALSE,square=FALSE)
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
