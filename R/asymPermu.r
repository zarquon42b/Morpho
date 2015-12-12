#' Assess differences in amount and direction of asymmetric variation (only object symmetry)
#'
#' Assess differences in amount and direction of asymmetric variation (only object symmetry)
#'
#' @param x object of class symproc result from calling \code{\link{procSym}} with \code{pairedLM} specified
#' @param groups factors determining grouping.
#' @param rounds number of permutations
#' @param which select which factorlevels to use, if NULL, all pairwise differences will be assessed after shuffling pooled data.
#' @note
#' This test is only sensible if between-group differences concerning directional asymmetry have been established (e.g. by applying a MANOVA on the "asymmetric" PCscores (see also \code{\link{procSym}}) and one wants to test whether these can be attributed to differences in amount and/or direction of asymmetric displacement. If there is no or only very little directional asymmetry present, the angles will only be significan when larger than 90 degrees (pi/2). So careful interpretation is advised. The Null-Hypothesis is that we have the same directional asymmetry in both groups. If you want to test whether the angle between groups is similar, please use \code{\link{angleTest}}.
#' @return
#' \item{dist }{difference between vector lengths of group means}
#' \item{angle }{angle (in radians) between vectors of group specific asymmetric deviation}
#' \item{means }{actual group averages}
#' \item{p.dist }{p-value obtained by comparing the actual distance to randomly acquired distances}
#' \item{p.angle }{p-value obtained by comparing the actual angle to randomly acquired angles}
#'  \item{permudist }{vector containing differences between random group means' vector lenghts}
#'  \item{permuangle }{vector containing angles between random group means' vectors}
#' \item{groupmeans}{ array with asymmetric displacement per group}
#' \item{levels}{ character vector containing the factors used}
#'
#' @seealso \code{\link{procSym}}
#' @export 
asymPermute <- function(x,groups,rounds=1000,which=NULL) {

    if (!inherits(x,"symproc"))
        stop("please provide object of class 'symproc'")
        
    if (is.list(x))
        asym <- vecx(x$Asym)
    else
        asym <- x
    groups <- factor(groups)
    class(asym) <- "symproc"
    lev <- levels(groups)
   
    if (!is.null(which)) {
        use <- which(groups %in% lev[which])
        groups <- factor(groups[use])
        lev <- levels(groups)
        asym <- asym[use,]
    }
    ng <- length(lev)
    if (length(groups) != nrow(asym))
        stop("number of groups and number of observations differ")
    gmeans <- matrix(NA, ng,ncol(asym))
    groupmeans <- array(NA, dim=c(dim(x$mshape),ng))
    for ( i in 1:ng) {
        gmeans[i,] <- colMeans(asym[groups == lev[i],,drop=FALSE])
        groupmeans[,,i] <- matrix(gmeans[i,],nrow(x$mshape), ncol(x$mshape))
    }
    
    shaker <- .Call("asymPermute",asym,as.integer(groups),as.integer(rounds))
    out <- list(groupmeans=groupmeans)
    out$levels <- lev
    if (rounds > 0) {
        dist <- matrix(0,ng,ng); dimnames(dist) <- list(lev,lev)
        angdiff <- dist.probs <- ang.probs <- dist
        count <- 1
        for (j1 in 1:(ng - 1)) {
            for (j2 in (j1 + 1):ng) {
                dist0 <- dist[j1,j2] <- shaker$dists[[count]][1]
                dists <- shaker$dists[[count]][-1]
                ang0 <- angdiff[j1,j2] <- shaker$angle[[count]][1]
                angs <- shaker$angles[[count]][-1]
                pdist.value <- length(which(dists >= dist0))
                pang.value <- length(which(angs >= ang0))
                if (pdist.value > 0) {
                    pdist.value <- pdist.value/rounds
                } else {
                    pdist.value <- 1/rounds
                }                
                dist.probs[j1,j2] <- pdist.value
                 if (pang.value > 0) {
                    pang.value <- pang.value/rounds
                } else {
                    pang.value <- 1/rounds
                }                
                dist.probs[j1,j2] <- pdist.value
                ang.probs[j1,j2] <- pang.value
                count <- count+1
            }
        }
        out$dist <- as.dist(t(dist)+dist)
        out$angle <- as.dist(angdiff+t(angdiff))
        out$p.dist <- as.dist(dist.probs+t(dist.probs))
        out$p.angle <- as.dist(ang.probs+t(ang.probs))
                                        #out$permudist <- sample.dists
                                        #out$permuangle <- sample.angle
    }
    
    return(out)
}
