#' @rdname asymPermute
#' @details asymPermuteArma is a version asymPermute making use of the RcppArmadillo package to speed up permutation.
#' @export asymPermuteArma
asymPermuteArma <- function(x,groups,rounds=1000,which=1:2) {

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
    mean1 <- meanMat(asym[groups == lev[1],])
    mean2 <- meanMat(asym[groups == lev[2],])

    shaker <- .Call("asymPerm",asym,as.integer(groups),as.integer(rounds))
    l.diff <- shaker$diff[1]
    a.diff <- shaker$angle[1]
    out <- list()
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
