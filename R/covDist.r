#' calculates distances and PC-coordinates of covariance matrices
#' 
#' calculates PC-coordinates of covariance matrices by using the Riemannian
#' metric in their respective space.
#' 
#' @details
#' \code{covDist} calculates the Distance between covariance matrices while \code{covPCA} uses a MDS (multidimensional scaling) approach to obtain PC-coordinates
#' from a  distance matrix derived from multiple groups. P-values for pairwise
#' distances can be computed by permuting group membership and comparing actual
#' distances to those obtained from random resampling. To calculate confidence intervals for PC-scores, within-group bootstrapping can be performed.
#' 
#' @param s1 m x m covariance matrix 
#' @param s2 m x m covariance matrix 
#' @return \code{covDist} returns the distance between s1 and s2
#' @author Stefan Schlager
#' @seealso \code{\link{prcomp}}
#' @references Mitteroecker P, Bookstein F. 2009. The ontogenetic trajectory of
#' the phenotypic covariance matrix, with examples from craniofacial shape in
#' rats and humans. Evolution 63:727-737.
#'
#' Hastie T, Tibshirani R, Friedman JJH.  2013. The elements of statistical
#' learning. Springer New York.
#' 
#'
#' @examples
#' 
#' 
#' cpca <- covPCA(iris[,1:4],iris[,5])
#' cpca$p.matrix #show pairwise p-values for equal covariance matrices
#' \dontrun{
#' require(car)
#' sp(cpca$PCscores[,1],cpca$PCscores[,2],groups=levels(iris[,5]),
#'    smooth=FALSE,xlim=range(cpca$PCscores),ylim=range(cpca$PCscores))
#' 
#' data(boneData)
#' proc <- procSym(boneLM)
#' pop <- name2factor(boneLM, which=3)
#' ## compare covariance matrices for PCscores of Procrustes fitted data
#' cpca1 <- covPCA(proc$PCscores, groups=pop, rounds = 1000)
#' ## view p-values:
#' cpca1$p.matrix # differences between covariance matrices
#' # are significant
#' ## visualize covariance ellipses of first 5 PCs of shape
#' spm(proc$PCscores[,1:5], groups=pop, smooth=FALSE,ellipse=TRUE, by.groups=TRUE)
#' ## covariance seems to differ between 1st and 5th PC
#' ## for demonstration purposes, try only first 4 PCs
#' cpca2 <- covPCA(proc$PCscores[,1:4], groups=pop, rounds = 1000)
#' ## view p-values:
#' cpca2$p.matrix # significance is gone
#' }
#'
#' @rdname covDist
#' @export
covDist <- function(s1,s2)
{
    dims1 <- dim(s1);dims2 <- dim(s2)
    
    if (dims1[1] != dims1[2] || dims2[1] != dims2[2] || dims1[1] != dims2[1])
        stop("please provide covariance matrices with idential dimensionality")
    cdist <- sqrt(sum(log(eigen(solve(s1,s2))$values)^2))
    return(cdist)
}

#' @param data matrix containing data with one row per observation
#' @param groups factor: group assignment for each specimen
#' @param rounds integer: rounds to run permutation of distances by randomly assigning group membership
#' @param bootrounds integer: perform bootstrapping to generate confidence intervals  (lower boundary, median and upper boundary) for PC-scores.
#' @param lower.bound numeric: set probability (quantile) for lower boundary estimate from bootstrapping.
#' @param upper.bound numeric: set probability (quantile) for upper boundary estimate from bootstrapping.
#' @return
#' 
#' \code{covPCA} returns a list containing:
#' 
#' if \code{scores = TRUE}
#' \item{PCscores }{PCscores}
#' \item{eigen}{eigen decomposition of the centered inner product}
#' if \code{rounds > 0}
#' \item{dist }{distance matrix}
#' \item{p.matrix }{p-values for pairwise distances from permutation testing}
#' if \code{bootrounds > 0}
#' \item{bootstrap}{list containing the lower and upper bound of the confidence intervals of PC-scores as well as the median of bootstrapped values.}
#' \item{boot.data}{array containing all results generated from bootstrapping.}
#'
#' @examples
#' #do some bootstrapping 1000 rounds
#' cpca <- covPCA(iris[,1:4],iris[,5],rounds=0, bootrounds=1000)
#' #plot bootstrapped data of PC1 and PC2 for first group
#' plot(t(cpca$boot.data[1,1:2,]),xlim=range(cpca$boot.data[,1,]),
#'                                ylim=range(cpca$boot.data[,2,]))
#' points(t(cpca$PCscores[1,]),col="white",pch=8,cex=1.5)##plot actual values
#'                       
#' for (i in 2:3) {
#'   points(t(cpca$boot.data[i,1:2,]),col=i)##plot other groups
#'   points(t(cpca$PCscores[i,]),col=1,pch=8,cex=1.5)##plot actual values
#' }
#' 
#' 
#' @rdname covDist
#' @export
#'
covPCA <- function(data,groups,rounds=1000,bootrounds=0,lower.bound=0.05, upper.bound=0.95) {
    out <- list()
    if (! is.factor(groups))
        groups <- as.factor(groups)
    
    groups <- factor(groups)
    lev <- levels(groups)
    nlev <- length(lev)
    for (i in 1:nlev) # center data per group
        data[groups==lev[i],] <- sweep(data[groups==lev[i],],2, colMeans(data[groups==lev[i],,drop=FALSE]))
    data <- as.matrix(data)
    groups <- as.integer(groups)
    rounds <- as.integer(rounds)[1]
    bootrounds <- as.integer(bootrounds)[1]
    result <- .Call("covPCAwrap", data, groups,bootrounds,rounds)
    
    dist <- result$dist
    dimnames(dist) <- list(lev,lev)
    out$dist <- as.dist(dist)
    ## permutation testing evaluation
    p.matrix <- matrix(NA, nlev, nlev)
    if (rounds > 0) {
        dist.mat <- result$permute
        for (i in 1:(nlev-1)) {
            for (j in (i+1):(nlev)) {
                sorti <- sort(dist.mat[j, i,  ])
                if (max(sorti) < dist[j, i]) {
                    p.matrix[j, i] <-  1/rounds
                } else {
                    marg <- min(which(sorti >= dist[j, i]))
                    p.matrix[j, i] <- (rounds - (marg-1))/rounds
                }
            }
        }
        dimnames(p.matrix) <- list(lev,lev)
        out$p.matrix <- as.dist(p.matrix)
    }
    
    out$PCscores <- result$Scores$PCscores
    out$eigen <- list(values=as.vector(result$Scores$eigenval), vectors=result$Scores$eigenvec)
    out$Var <- out$eigen$values/sum(out$eigen$values)
    if (bootrounds > 0) {
        bootres <- result$bootstrap
        quants05 <- apply(bootres,1:2,quantile,probs=lower.bound)
        quants95 <- apply(bootres,1:2,quantile,probs=upper.bound)
        quants50 <- apply(bootres,1:2,quantile,probs=0.5)
        boot <- list()
        boot$lower.bound <- quants05
        boot$upper.bound <- quants95
        boot$median <- quants50
        out$boot.data <- bootres
        out$bootstrap <- boot
    }
        
    
    return(out)
}
