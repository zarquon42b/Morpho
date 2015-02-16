.CVAdists <- function(N, groups,rounds,winv,p.adjust.method) {
    if (rounds == 0)
        rounds <- 1
    l <- dim(N)[2]
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
    
    groups <- as.integer(groups)
    shaker <- .Call("CVAdists", N, groups, as.integer(rounds), winv)
    mahadist <- mahaprobs <- procdist <- procprobs <- matrix(0,ng,ng)
    if (!is.null(lev)) 
        dimnames(mahadist) <- dimnames(mahaprobs) <- dimnames(procdist) <- dimnames(procprobs) <- list(lev,lev)
    count <- 1
    for (j1 in 1:(ng - 1)) {
        for (j2 in (j1 + 1):ng) {
            ## check plain distances
            dist0 <- procdist[j1,j2] <- shaker$Plain[[count]][1]
            dists <- shaker$Plain[[count]][-1]
            p.value <- length(which(dists >= dist0))
            if (p.value > 0) {
                p.value <- p.value/rounds
            } else {
                p.value <- 1/rounds
            }
            procprobs[j1,j2] <- p.value

            dist0 <- mahadist[j1,j2] <- shaker$Maha[[count]][1]
            dists <- shaker$Maha[[count]][-1]
            p.value <- length(which(dists >= dist0))
            if (p.value > 0) {
                p.value <- p.value/rounds
            } else {
                p.value <- 1/rounds
            }
            mahaprobs[j1,j2] <- p.value
            
            count <- count+1
        }
    }
    out <- list()
    out$GroupdistMaha <- as.dist(t(mahadist))
    out$GroupdistEuclid <- as.dist(t(procdist))
    if (rounds > 1) {
        mahaprobs[upper.tri(mahaprobs)] <- p.adjust(mahaprobs[upper.tri(mahaprobs)],method=p.adjust.method)
        procprobs[upper.tri(procprobs)] <- p.adjust(procprobs[upper.tri(procprobs)],method=p.adjust.method)
        out$probsMaha <- as.dist(t(mahaprobs))
        out$probsEuclid <- as.dist(t(procprobs))
        out$p.adjust.method <- match.arg(p.adjust.method,p.adjust.methods)
    }
    return(out)
}



