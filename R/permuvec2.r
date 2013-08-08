permuvec2 <- function(data, groups, mc.cores=detectCores(), rounds=0)
{
    n3 <- FALSE
    if (length(dim(data))==3) {
        data <- vecx(data)
        n3 <- TRUE
    }
    groups <- factor(groups)
    lev <- levels(groups)
    nlev <- length(lev)
    Gmeans <- apply(data[groups==lev[[1]],],2,mean)
    for (i in 2:nlev)
        Gmeans <- rbind(Gmeans, apply(data[groups==lev[[i]],],2,mean))

    vlen <- apply(Gmeans, 1, function(x) x <- sqrt(sum(x^2)))
    vecdiff <- as.matrix(dist(vlen))
    anglemat <- matrix(0, nlev, nlev)
    for (i in 1:(nlev - 1)) {
        for (j in (i + 1):nlev) 
            anglemat[j, i] <- angle.calc(Gmeans[i, ], Gmeans[j,])
    }
    rownames(anglemat) <- rownames(vecdiff) <- lev
    out <- list()
    out$angle <- as.dist(anglemat)
    out$vecdiff <- as.dist(vecdiff)
    if (rounds > 0) {
        registerDoParallel(mc.cores)
        p.vecdiff <- p.angle <- matrix(0, nlev, nlev)
        rownames(p.vecdiff) <- rownames(p.angle) <- lev
        p.vecdiff.arr <- p.angle.arr <- array(NA, dim=c(nlev, nlev, rounds))
        permudist <- function(i)
            {
                samplegroup <- sample(groups)
                out <- permuvec2(data, samplegroup)
                return(out)
            }
        if(.Platform$OS.type == "windows")
            outlist <- foreach (i = 1:rounds) %do% permudist(i)
        else
            outlist <- foreach (i = 1:rounds) %dopar% permudist(i)
        for (i in 1:rounds) {
            p.vecdiff.arr[,,i] <- as.matrix(outlist[[i]]$vecdiff)
            p.angle.arr[,,i] <- as.matrix(outlist[[i]]$angle)
        }
        for (i in 1:(nlev - 1)) {
            for (j in (i + 1):nlev) {
                sortv <- sort(p.vecdiff.arr[j, i, ])
                if (max(sortv) < as.matrix(vecdiff[j, i])) {
                    p.vecdiff[j, i] <- 1/rounds
                } else {
                    marg <- min(which(sortv >= vecdiff[j, i]))
                    p.vecdiff[j, i] <- (rounds - marg+1)/rounds
                }
                sorta <- sort(p.angle.arr[j, i, ])
                if (max(sorta) < as.matrix(anglemat[j, i])) {
                    p.angle[j, i] <- 1/rounds
                } else {
                    marg <- min(which(sorta >= anglemat[j, i]))
                    p.angle[j, i] <- (rounds - marg+1)/rounds
                }
            }
        }
        out$p.vecdiff <- as.dist(p.vecdiff)
        out$p.angle <- as.dist(p.angle)
    }
    return(out)
}



