.CVAcrova <- function(dataarray,groups,test ,weighting=TRUE,tolinv=1e-10,robust=c("classical", "mve", "mcd"),...)
{   
    groups <- factor(groups)
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
   
    N <- dataarray
    n3 <- FALSE
    N <- as.matrix(N)
    n <- dim(N)[1]
    l <- dim(N)[2]
    if (length(groups) != n)
        warning("group affinity and sample size not corresponding!")
    
    if (is.vector(N) || dim(N)[2] == 1)
        stop("data should contain at least 2 variable dimensions")
    
    covWithin <- covW(N, groups,robust=robust,...)
    Gmeans <- attributes(covWithin)$means
    if (weighting) {
        Grandm <- colSums(Gmeans*gsizes)/n ## calculate weighted Grandmean (thanks to Anne-Beatrice Dufour for the bug-fix)
    } else {
        Grandm <- colMeans(Gmeans)
    }
    N <- sweep(N, 2, Grandm) #center data according to Grandmean
    resGmeans <- sweep(Gmeans, 2, Grandm)
    
    if (weighting) {
        for (i in 1:ng) 
            resGmeans[i, ] <- sqrt(gsizes[i]) * resGmeans[i, ]
        X <- resGmeans
    } else 
        X <- sqrt(n/ng) * resGmeans
    
    eigcoW <- eigen(covWithin); ## eigen decomp of between group covariance Matrix
    E <- eigcoW$values*(n - ng)  ##eigenvalues of between group SSPQR
    U <- eigcoW$vectors
    Ec <- eigcoW$values
    Ec2 <- Ec
    geninv <- FALSE
    if (min(Ec) < tolinv) {
        #cat(paste("singular Covariance matrix: General inverse is used. Threshold for zero eigenvalue is", tolinv, "\n"))
        geninv <- TRUE
    }
    abovetol <- which(Ec >= tolinv)
    E[abovetol] <- sqrt(1/E[abovetol])
    Ec[abovetol] <- sqrt(1/Ec[abovetol])
    Ec2[abovetol] <- 1/Ec2[abovetol]
    if (geninv)
        Ec[-abovetol] <- E[-abovetol] <- Ec2[-abovetol] <- 0

    useEig <- min((ng-1), l)
    ZtZ <- (E * t(X %*% U))
    eigZ <- svd(ZtZ,nv=0,nu=useEig)
    eigZ$d <- eigZ$d^2
    A <- Re(eigZ$u)
    CV <- U %*% (Ec * A)
    di <- dim(CV)[2]
    
    for (i in 1:di) {
        rho <- angle.calc(test[,i],CV[,i])
        if (rho > pi/2)
            CV[,i] <- -CV[,i]		
    }
    meanscores <-  sweep(Gmeans, 2, Grandm) %*%CV
    rownames(meanscores) <- lev
    return(list(CV=CV,Grandmean=Grandm,meanscores=meanscores))
}
