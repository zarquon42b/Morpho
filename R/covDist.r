covDist <- function(s1,s2)
{
    dims1 <- dim(s1);dims2 <- dim(s2)
    
    if (dims1[1] != dims1[2] || dims2[1] != dims2[2] || dims1[1] != dims2[1])
        stop("please provide covariance matrices with idential dimensionality")
    cdist <- sqrt(sum(log(eigen(solve(s1,s2))$values)^2))
    return(cdist)
}
covPCA <- function(data,groups,scores=TRUE,rounds=0, mc.cores=detectCores())
{
    out <- list()
    if (! is.factor(groups))
        groups <- as.factor(groups)
    
    groups <- factor(groups)
    lev <- levels(groups)
    nlev <- length(lev)
    for (i in 1:nlev) # center data per group
        data[groups==lev[i],] <- sweep(data[groups==lev[i],],2, apply(data[groups==lev[i],],2,mean))
    covlist <- list()
    for (i in 1:nlev)
        covlist[[i]] <- cov(data[groups==lev[i],])
    V <- diag(0,nlev,nlev)
    for (i in 1:(nlev-1)) {
        for (j in (i+1):(nlev))
            V[j,i] <- covDist(covlist[[j]],covlist[[i]])^2
    }
    V <- V+t(V)
    dimnames(V) <- list(lev,lev)
    out$dist <- as.dist(V)
    if (rounds > 0)
        out$p.matrix <- .covPCApermut(data, groups, rounds, mc.cores, V)
    
    if (scores) {
        H <- matrix(-1/nlev,nlev,nlev)
        H <- H+diag(nlev)
        D <- (-1/2)*(H%*%V%*%H)
        eigenD <- eigen(D,symmetric = TRUE)
        eigenD$values <- eigenD$values[1:(nlev-1)]
        eigenD$vectors <- eigenD$vectors[,1:(nlev-1)]
        PCscores <- as.matrix(t(t(eigenD$vectors)*sqrt(eigenD$values)))
        rownames(PCscores) <- lev
        
        out$PCscores <- PCscores
        out$Var <- eigenD$values/sum(eigenD$values)
        out$eigen <- eigenD
    }
    return(out)
}

.covPCApermut <- function(data, groups, rounds, mc.cores, V)
{
    if(.Platform$OS.type == "windows")
        mc.cores=1
    lev <- levels(groups)
    nlev <- length(levels(groups))
    p.matrix <- matrix(NA, nlev, nlev)
    dist.mat <- array(0, dim = c(nlev, nlev, rounds))
    registerDoParallel(cores=mc.cores)
    permufun <- function(i)
        {
            permugroup <- sample(groups)
            distout <- as.matrix(covPCA(data, permugroup, scores = FALSE)$dist)
            return(distout)
        }
    if(.Platform$OS.type == "windows")
        a.list <- foreach(i=1:rounds)%do%permufun(i)
    else
        a.list <- foreach(i=1:rounds)%dopar%permufun(i)
    for (i in 1:rounds)
        dist.mat[,,i] <- a.list[[i]]
    for (i in 1:(nlev-1)) {
        for (j in (i+1):(nlev)) {
            sorti <- sort(dist.mat[j, i,  ])
            if (max(sorti) < V[j, i]) {
                p.matrix[j, i] <-  1/rounds
            } else {
                marg <- min(which(sorti >= V[j, i]))
                p.matrix[j, i] <- (rounds - (marg-1))/rounds
            }
        }
    }
    dimnames(p.matrix) <- list(lev,lev)
    return(as.dist(p.matrix))
}

