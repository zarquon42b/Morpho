.CVAdists <- function(N, Tmatrix, lev, rounds, proc.dist, maha.dist, mc.cores,n3, winv,groups )
{

    n <- dim(N)[1]
    l <- dim(N)[2]
    ng <- length(groups)
    pmatrix.proc <- pmatrix <- matrix(NA, ng, ng) ### generate distance matrix Mahal
    
    if (!is.null(lev))
        {
            rownames(pmatrix) <- lev
            colnames(pmatrix) <- lev
            rownames(pmatrix.proc) <- lev
            colnames(pmatrix.proc) <- lev
        }
    if (rounds > 0)
        {
            roun <- function(i)
                {   
                    groups1 <- list(numeric(0))
                    dist.mat <- matrix(0,ng,ng)
                    shake <- sample(1:n)
                    Gmeans1 <- matrix(0, ng, l)
                    l1 <- 0
                    
                    for (j in 1:ng)
                        {
                            groups1[[j]] <- c(shake[(l1 + 1):(l1 + (length(groups[[j]])))])
                            l1 <- l1 + length(groups[[j]])
                            tmpmat <- N[groups1[[j]],]
                            if ( length(groups[[j]])==1)
                                tmpmat <- t(as.matrix(tmpmat))
                            Gmeans1[j, ] <- apply(tmpmat, 2, mean)
                        }
                    for (j1 in 1:(ng - 1)) 
                        for (j2 in (j1 + 1):ng) 
                            dist.mat[j2, j1] <- sqrt((Gmeans1[j1, ] - Gmeans1[j2, ]) %*% winv %*% (Gmeans1[j1,] - Gmeans1[j2, ]))
                    return(dist.mat)
                }
            
            dist.mat <- array(0, dim = c(ng, ng, rounds))
            a.list <- as.list(1:rounds)
            a.list <- mclapply(a.list,roun,mc.cores=mc.cores)
            
            for (i in 1:rounds)
                dist.mat[,,i] <- a.list[[i]]
            
            for (j1 in 1:(ng - 1)) 
                {
                    for (j2 in (j1 + 1):ng) 
                        {
                            sorti <- sort(dist.mat[j2, j1, ])
                            if (max(sorti) < maha.dist[j2, j1]) 
                                pmatrix[j2, j1] <- 1/rounds
                            else 
                                {
                                    marg <- min(which(sorti >= maha.dist[j2, j1]))
                                    pmatrix[j2, j1] <- (rounds - marg+1)/rounds
                                }
                        }
                }
            
            roun.proc <- function(i)
                {
                    groups1 <- list()
                    dist.mat <- matrix(0,ng,ng)
                    shake <- sample(1:n)
                    Gmeans1 <- matrix(0, ng, l)
                    l1 <- 0
                    for (j in 1:ng) 
                        {
                            groups1[[j]] <- c(shake[(l1 + 1):(l1 + (length(groups[[j]])))])
                            l1 <- l1 + length(groups[[j]])
                            tmpmat <- Tmatrix[groups1[[j]],]
                            if ( length(groups[[j]]) > 1)
                                Gmeans1[j, ] <- apply(tmpmat, 2, mean)
                            else 
                                Gmeans1[j, ] <- tmpmat
                        }
                    
                    for (j1 in 1:(ng - 1)) 
                        for (j2 in (j1 + 1):ng)
                            {
                                if (n3)
                                    dist.mat[j2, j1] <- angle.calc(Gmeans1[j1, ],Gmeans1[j2, ])
                                else
                                    dist.mat[j2, j1] <- sqrt(sum((Gmeans1[j1, ]-Gmeans1[j2, ])^2))   
                                return(dist.mat)
                            }
                }
            
            dist.mat.proc <- array(0, dim = c(ng, ng, rounds))
            a.list <- as.list(1:rounds)
            a.list <- mclapply(a.list,roun.proc,mc.cores=mc.cores)
            
            for (i in 1:rounds)
                dist.mat.proc[,,i] <- a.list[[i]]
            
            for (j1 in 1:(ng - 1)) 
                {
                    for (j2 in (j1 + 1):ng) 
                        {
                            sorti <- sort(dist.mat.proc[j2, j1, ])
                            if (max(sorti) < proc.dist[j2, j1]) 
                                pmatrix.proc[j2, j1] <- 1/rounds
                            else 
                                {
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
