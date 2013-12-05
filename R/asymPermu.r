#' Assess differences in amount and direction of asymmetric variation
#'
#' Assess differences in amount and direction of asymmetric variation
#'
#' @param x object of class symproc result from calling \code{\link{procSym}} with \code{pairedLM} specified
#' @param groups factors determining grouping.
#' @param rounds number of permutations
#' @param which in case the factor levels are >2 this determins which
#' factorlevels to use
#' @param mc.cores integer: determines how many cores to use for the
#' computation. The default is autodetect. But in case, it doesn't work as
#' expected cores can be set manually. Parallel processing is disabled on
#' Windows due to occasional errors.
#' @return
#' \item{dist }{difference between vector lengths of group means}
#' \item{angle }{angle between vectors of group specific asymmetric deviation}
#' \item{means }{actual group averages}
#' \item{p.dist }{p-value obtained by comparing the actual distance to randomly acquired distances}
#' \item{p.angle }{p-value obtained by comparing the actual angle to randomly acquired angles}
#'  \item{permudist }{vector containing differences between random group means' vector lenghts}
#'  \item{permuangle }{vector containing angles between random group means' vectors}
#'
#' @importFrom foreach registerDoSEQ
#' @export asymPerm
#' 
asymPerm <- function(x,groups,rounds=1000,which=1:2,mc.cores=detectCores()) {

    if (!inherits(x,"symproc"))
        stop("please provide object of class 'symproc'")
    if (mc.cores > 1 && rounds > 0) {
        if (.Platform$OS.type != "windows" ) {
            registerDoParallel(cores=mc.cores)
        } else {
            cl <- makeCluster(mc.cores)
            registerDoParallel(cl=cl)
        }
    } else if (mc.cores == 1 && rounds > 0)
        registerDoSEQ()
    
    if (is.list(x))
        asym <- vecx(x$Asym)
    else
        asym <- x
    class(asym) <- "symproc"
    lev <- levels(groups)
    if (length(lev) > 2) {
        use <- which(groups %in% lev[which])
        groups <- factor(groups[use])
        lev <- levels(groups)
        asym <- asym[use,]
    }
    mean1 <- meanMat(asym[groups == lev[1],])
    mean2 <- meanMat(asym[groups == lev[2],])
    l.diff <- abs(sqrt(sum(mean1^2))-sqrt(sum(mean2^2)))
    a.diff <- angle.calc(mean1,mean2)
    out <- list()
    out$dist <- l.diff
    out$angle <- a.diff
    permuta <- function() {
        shake <- sample(groups)
        out <- asymPerm(asym,shake,rounds=0,which=which,mc.cores=1)
        return(out)
    }
    if (rounds > 0) {
        i=0
        resampled <- foreach(i = 1:rounds,.inorder = FALSE,.packages=c("Morpho")) %dopar% permuta()
        out$means <- rbind(mean1,mean2)
        rownames(out$means) <- lev
        even <- 1:(rounds)*2
        odd <- even-1
        resampled <- unlist(resampled)
        sample.dists <- resampled[odd]
        sample.angle <- resampled[even]
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
        if (.Platform$OS.type == "windows" && mc.cores > 1)
            stopCluster(cl)
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
#' @export meanMat
meanMat <-function(A,usedim=2)
    {
        if (usedim==1)
            A <- t(A)
        out <- A[1,]*0
        resul <- .Fortran("meanMat",A,nrow(A),ncol(A),out=out)
        out <- resul$out
        return(out)
    }
