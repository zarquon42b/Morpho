.groupPCAcrova <- function(N, groups,tol=1e-10,groupPCs,weighting=weighting)
{ 
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
    
    if (length(dim(N)) == 3) 
        N <- vecx(N)
    #n <- dim(N)[1]
    l <- dim(N)[2]
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng) {
        Gmeans[i, ] <- colMeans(N[groups==lev[i], ,drop=F])
    }
    if (weighting)
        wt <- gsizes
    else
    wt <- rep(1,ng)
    wcov <- cov.wt(Gmeans,wt=wt)
    Grandm <- as.vector(wcov$center)
    eigenGmeans <- eigen(wcov$cov)
    N <- sweep(N, 2, Grandm)
    valScores <- which(eigenGmeans$values > tol)
    #groupScores <- N%*%(eigenGmeans$vectors[,valScores])
    PCs <- as.matrix(eigenGmeans$vectors[,valScores])
    
        if (is.null(dim(groupPCs))) {
            PCs <- matrix(PCs,length(PCs),1)
            groupPCs <- matrix(groupPCs,length(groupPCs),1)
        }
    di <- dim(PCs)[2]
    for (i in 1:di) {
        rho <- angle.calc(groupPCs[,i ],PCs[,i])
        if (rho > pi/2)
            PCs[,i] <- PCs[,i]*(-1)
    }
    return(list(PCs=PCs,Grandmean=Grandm))
}
