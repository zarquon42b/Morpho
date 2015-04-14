#' determine the minimum ratio for two subsequent eigenvalues to be considered different
#'
#' determine the minimum ratio for two subsequent eigenvalues to be considered different
#' @param n sample size
#' @param expect expectation value for chi-square distribution of df=2
#' @return returns the minimum ratio between two subsequent subsequent eigenvalues to be considered different.
#' @examples
#' ## reproduce the graph from Bookstein (2014, p. 324)
#' ## and then compare it to ratios for values to be considered
#' ## statistically significant
#' myseq <- seq(from=10,to = 50, by = 2)
#' myseq <- c(myseq,seq(from=50,to=1000, by =20))
#' ratios <- getPCtol(myseq)
#' plot(log(myseq),ratios,cex=0,xaxt = "n",ylim=c(1,5.2))
#' ticks <- c(10,20,50,100,200,300,400,500,600,700,800,900,1000)
#' axis(1,at=log(ticks),labels=ticks)
#' lines(log(myseq),ratios)
#' abline(v=log(ticks), col="lightgray", lty="dotted")
#' abline(h=seq(from=1.2,to=5, by = 0.2), col="lightgray", lty="dotted")
#'
#' ## now we raise the bar and compute the ratios for values
#' ## to be beyond the 95th percentile of
#' ## the corresponding chi-square distribution:
#' ratiosSig <- getPCtol(myseq,expect=qchisq(0.95,df=2))
#' lines(log(myseq),ratiosSig,col=2)
#'
#' 
#' @references
#' Bookstein, F. L. Measuring and reasoning: numerical inference in the sciences. Cambridge University Press, 2014
#' @seealso  \code{\link{getMeaningfulPCs}}
#' @export
getPCtol <- function(n,expect=2) {
    x <- exp(expect/(2*n))+ sqrt(exp((2*expect)/(2*n))-1)
    return(x^2)
}
#' get number of meaningful Principal components
#'
#' get number of meaningful Principal components
#' @param values eigenvalues from a PCA
#' @param n sample size
#' @param expect expectation value for chi-square distribution of df=2
#' @param sdev logical: if TRUE, it is assumed that the values are square roots of eigenvalues.
#' @return
#' \item{tol}{threshold of ratio specific for \code{n}}
#' \item{good}{integer vector specifying the meaningful Principal Components}
#' @details This implements the method suggested by Bookstein (2014, pp. 324), to determine whether a PC is entitled to interpretation. I.e. a PC is regarded meaningful (its direction) if the ratio of this PC and its successor is above a threshold based on a log-likelihood ratio (and dependend on sample size).
#' @examples
#' data(boneData)
#' proc <- procSym(boneLM)
#' getMeaningfulPCs(proc$eigenvalues,n=nrow(proc$PCscores))
#' ## the first 3 PCs are reported as meaningful
#' ## show barplot that seem to fit the bill
#' barplot(proc$eigenvalues)
#' @references
#' Bookstein, F. L. Measuring and reasoning: numerical inference in the sciences. Cambridge University Press, 2014
#' @seealso \code{\link{getPCtol}}
#' @export
getMeaningfulPCs <- function(values,n,expect=2,sdev=FALSE) {
    nv <- length(values)
    if (sdev)
        values <- values^2
    tol <- getPCtol(n,expect)
    coeff <- values[1:(nv-1)]/values[2:nv]
    
    bad <- which(coeff < tol)
    out <- list(tol=tol)
    mbad <- min(bad)
    if (mbad < 2) {
        warning("there are no meaningful PCs")
        out$good <- NULL
    } else
        out$good <- (1:nv)[1:(mbad-1)]
    return(out)
}
    
