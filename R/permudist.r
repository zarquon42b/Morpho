#' performs permutation testing for group differences.
#' 
#' This function compares the distance between two groupmeans to the distances
#' obtained by random assignment of observations to this groups.
#' 
#' 
#' @param data array or matrix containing data
#' @param groups factors determining grouping.
#' @param rounds number of permutations
#' @param which integer (optional): in case the factor levels are > 2 this determins which factorlevels to use
#' @param p.adjust.method method to adjust p-values for multiple comparisons see \code{\link{p.adjust.methods}} for options.
#' 
#' @return
#' \item{dist }{distance matrix with distances between actual group means}
#' \item{p.adjust.method}{method used for p-value adjustion}
#' \item{p.value }{distance matrix containing pairwise p-values obtained by comparing the actual distance to randomly acquired distances}
#' 
#' 
#' @examples
#' 
#' data(boneData)
#' proc <- procSym(boneLM)
#' groups <- name2factor(boneLM,which=3)
#' perm <- permudist(proc$PCscores[,1:10], groups=groups, rounds=10000)
#' 
#' ## now we concentrate only on sex dimorphism between Europeans
#' groups <- name2factor(boneLM,which=3:4)
#' levels(groups)
#' perm1 <- permudist(proc$PCscores, groups=groups,which=3:4, rounds=10000)
#' 
#' 
#' @export 
permudist <- function(data, groups, rounds=1000, which=NULL,p.adjust.method= "none")
{
    if (rounds == 0)
        rounds <- 1
    out <- list()
### configure grouping ####
    N <- data
    if (is.vector(N)) {
        N <- as.matrix(N)
    } else if (length(dim(N)) == 3)
        N <- vecx(N)
    
    if (!is.factor(groups))
        groups <- factor(groups)
    
    if (is.factor(groups)) {
        groups <- factor(groups)
        lev <- levels(groups)
                                        #levn <- length(lev)
    }
    old_groups <- groups
    if (!is.null(which)) {
        groups <- factor(groups[groups %in% lev[which]])
        lev <- levels(groups)
        N <- N[which(old_groups %in% lev),]
    }
    ng <- length(lev)
    N <- as.matrix(N)
    if (ng < 2)
        stop("provide at least two groups")
    if (length(groups) != nrow(N))
        warning("group affinity and sample size not corresponding!")
### end configure grouping ####
    if (rounds > 0) {
        shaker <- .Call("permudistArma",N,as.integer(groups),as.integer(rounds))
        dist <- matrix(0,ng,ng); dimnames(dist) <- list(lev,lev)
        probs <- dist
        count <- 1
        for (j1 in 1:(ng - 1)) {
            for (j2 in (j1 + 1):ng) {
                dist0 <- dist[j1,j2] <- shaker[[count]][1]
                dists <- shaker[[count]][-1]
                p.value <- length(which(dists >= dist0))
                if (p.value > 0) {
                    p.value <- p.value/rounds
                } else {
                    p.value <- 1/rounds
                }
                probs[j1,j2] <- p.value
                count <- count+1
            }
        }
        probs[upper.tri(probs)] <- p.adjust(probs[upper.tri(probs)],method=p.adjust.method)
        
        probs <- as.dist(probs+t(probs))
        dist <- as.dist(dist+t(dist))
        out$p.value <-probs
        out$p.adjust.method <- match.arg(p.adjust.method,p.adjust.methods)
    }
    out$dist <- dist

    return(out)
}
