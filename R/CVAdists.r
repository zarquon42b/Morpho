.CVAdists <- function(N, Tmatrix, groups, rounds, proc.dist, maha.dist, mc.cores,n3, winv )
{   ## calculates and permutes mahalanobis and euclidean distances for CVA
    #n <- dim(N)[1]
    l <- dim(N)[2]
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
    pmatrix.proc <- pmatrix <- matrix(NA, ng, ng) ### generate distance matrix Mahal
    
    if (!is.null(lev)) {
            rownames(pmatrix) <- lev
            colnames(pmatrix) <- lev
            rownames(pmatrix.proc) <- lev
            colnames(pmatrix.proc) <- lev
        }
    if (rounds > 0) {
        permuMaha <- function(i)
            {   
                dist.mat <- matrix(0,ng,ng)
                shake <- sample(groups)
                Gmeans.tmp <- matrix(0, ng, l)
                for (j in 1:ng) {
                    if(gsizes[j] > 1)
                        Gmeans.tmp[j, ] <- apply(N[shake==lev[j], ], 2, mean)
                    else
                        Gmeans.tmp[j, ] <- N[shake==lev[j], ]
                }
                for (j1 in 1:(ng - 1)) 
                    for (j2 in (j1 + 1):ng) 
                        dist.mat[j2, j1] <- sqrt((Gmeans.tmp[j1, ] - Gmeans.tmp[j2, ]) %*% winv %*% (Gmeans.tmp[j1,] - Gmeans.tmp[j2, ]))
                return(dist.mat)
            }
        
        dist.mat <- array(0, dim = c(ng, ng, rounds))
        a.list <- as.list(1:rounds)
        a.list <- mclapply(a.list, permuMaha, mc.cores=mc.cores)
        for (i in 1:rounds)
            dist.mat[,,i] <- a.list[[i]]
        
        for (j1 in 1:(ng - 1)) {
            for (j2 in (j1 + 1):ng) {
                sorti <- sort(dist.mat[j2, j1, ])
                if (max(sorti) < maha.dist[j2, j1]) 
                    pmatrix[j2, j1] <- 1/rounds
                else {
                    marg <- min(which(sorti >= maha.dist[j2, j1]))
                    pmatrix[j2, j1] <- (rounds - marg+1)/rounds
                }
            }
        }
        permuProc <- function(i)
            {
                dist.mat <- matrix(0,ng,ng)
                shake <- sample(groups)
                Gmeans.tmp <- matrix(0, ng, l)
                for (j in 1:ng) {
                    if(gsizes[j] > 1)
                        Gmeans.tmp[j, ] <- apply(N[shake==lev[j], ], 2, mean)
                    else
                        Gmeans.tmp[j, ] <- N[shake==lev[j], ]
                }
                for (j1 in 1:(ng - 1)) 
                    for (j2 in (j1 + 1):ng) {
                        if (n3)
                            dist.mat[j2, j1] <- angle.calc(Gmeans.tmp[j1, ],Gmeans.tmp[j2, ])
                        else
                            dist.mat[j2, j1] <- sqrt(sum((Gmeans.tmp[j1, ]-Gmeans.tmp[j2, ])^2))   
                        return(dist.mat)
                    }
            }
        dist.mat.proc <- array(0, dim = c(ng, ng, rounds))
        a.list <- as.list(1:rounds)
        a.list <- mclapply(a.list,permuProc,mc.cores=mc.cores)
        
        for (i in 1:rounds)
            dist.mat.proc[,,i] <- a.list[[i]]
        
        for (j1 in 1:(ng - 1)) {
            for (j2 in (j1 + 1):ng) {
                sorti <- sort(dist.mat.proc[j2, j1, ])
                if (max(sorti) < proc.dist[j2, j1]) 
                    pmatrix.proc[j2, j1] <- 1/rounds
                else {
                    marg <- min(which(sorti >= proc.dist[j2, j1]))
                    pmatrix.proc[j2, j1] <- (rounds - marg+1)/rounds
                }
            }
        }
    }
    pmatrix.proc <- as.dist(pmatrix.proc)
    pmatrix <- as.dist(pmatrix)
    maha.dist <- as.dist(maha.dist)
    proc.dist <- as.dist(proc.dist)
    Dist <- list(GroupdistMaha = maha.dist,GroupdistProc=proc.dist, probsMaha = pmatrix,probsProc = pmatrix.proc)
    return(Dist)
}
