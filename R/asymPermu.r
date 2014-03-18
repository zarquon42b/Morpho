#' Assess differences in amount and direction of asymmetric variation
#'
#' Assess differences in amount and direction of asymmetric variation
#'
#' @param x object of class symproc result from calling \code{\link{procSym}} with \code{pairedLM} specified
#' @param groups factors determining grouping.
#' @param rounds number of permutations
#' @param which in case the factor levels are >2 this determins which
#' factorlevels to use
#' @return
#' \item{dist }{difference between vector lengths of group means}
#' \item{angle }{angle between vectors of group specific asymmetric deviation}
#' \item{means }{actual group averages}
#' \item{p.dist }{p-value obtained by comparing the actual distance to randomly acquired distances}
#' \item{p.angle }{p-value obtained by comparing the actual angle to randomly acquired angles}
#'  \item{permudist }{vector containing differences between random group means' vector lenghts}
#'  \item{permuangle }{vector containing angles between random group means' vectors}
#' \item{groupmeans}{ array with asymmetric displacement per group}
#'
#'
#' @export 
asymPermute <- function(x,groups,rounds=1000,which=1:2) {

    if (!inherits(x,"symproc"))
        stop("please provide object of class 'symproc'")
        
    if (is.list(x))
        asym <- vecx(x$Asym)
    else
        asym <- x
    groups <- factor(groups)
    class(asym) <- "symproc"
    lev <- levels(groups)
    if (length(lev) > 2) {
        use <- which(groups %in% lev[which])
        groups <- factor(groups[use])
        lev <- levels(groups)
        asym <- asym[use,]
    }
    if (length(groups) != nrow(asym))
        stop("number of groups and number of observations differ")
    mean1 <- meanMat(asym[groups == lev[1],])
    mean2 <- meanMat(asym[groups == lev[2],])
    mean1Mat <- matrix(mean1,nrow(x$mshape), ncol(x$mshape))
    mean2Mat <- matrix(mean2,nrow(x$mshape), ncol(x$mshape))
    
    shaker <- .Call("asymPerm",asym,as.integer(groups),as.integer(rounds))
    l.diff <- shaker$diff[1]
    a.diff <- shaker$angle[1]
    out <- list()
    out$groupmeans <- bindArr(mean1Mat,mean2Mat,along=3)
    out$dist <- l.diff
    out$angle <- a.diff
    if (rounds > 0) {
        sample.dists <- shaker$diff[-1]
        sample.angle <- shaker$angle[-1]
        ## calculate p-values
        p.angle <- length(which(sample.angle >= a.diff))
        if (!!p.angle) {
            p.angle <- p.angle/rounds
            names(p.angle) <- "p-value"
        } else {
            p.angle <- 1/rounds
            names(p.angle) <- "p-value <"
        }
        p.dist <- length(which(sample.dists >= l.diff))
        if (!!p.dist) {
            p.dist <- p.dist/rounds
            names(p.dist) <- "p-value"
        } else {
            p.dist <- 1/rounds
            names(p.dist) <- "p-value <"
        }
        out$p.dist <- p.dist
        out$p.angle <- p.angle
        out$permudist <- sample.dists
        out$permuangle <- sample.angle
    }
        
    return(out)
}
#' fast calculation of a Matrix' per row/ per column mean - useful for very large matrices
#'
#' fast calculation of a Matrix' per row/ per column mean - equivalent to apply(X,2,mean) or apply(X,1,mean)- useful for very large matrices
#'
#' @param A numeric matrix
#' @param usedim integer: select over which dimension to average
#'
#' @return
#' vector containing row/column mean
#' @examples
#' A <- matrix(rnorm(1e6),1000,1000)
#' b <- meanMat(A)
#' # same as apply(A,2,mean)
#' b1 <- meanMat(A,1)
#' # same as apply(A,1,mean)
#' \dontrun{
#' #compare timing
#' system.time(meanMat(A))
#' system.time(apply(A,2,mean))
#' }
#' @export
meanMat <-function(A,usedim=2)
{
    if (!is.matrix(A) || !is.numeric(A))
        stop("A must be a numeric matrix")
    if (usedim==1)
        A <- t(A)
    out <- A[1,]*0
    resul <- .Call("meanMat",A)
    return(resul)
} 
