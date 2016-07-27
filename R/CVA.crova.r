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
    Tmatrix <- N
    N <- sweep(N, 2, Grandm) #center data according to Grandmean
    resGmeans <- sweep(Gmeans, 2, Grandm)
    
    if (weighting) {
        for (i in 1:ng) 
            resGmeans[i, ] <- sqrt(gsizes[i]) * resGmeans[i, ]
        X <- resGmeans
    } else 
        X <- sqrt(n/ng) * resGmeans
    
    covWithin <- covW(N, groups)
    eigW <- eigen(covWithin*(n - ng))
    eigcoW <- eigen(covWithin)
    U <- eigW$vectors
    E <- eigW$values
    Ec <- eigcoW$values
    Ec2 <- Ec
    
    #if (min(E) < tolinv)
    #    cat(paste("singular Covariance matrix: General inverse is used. Threshold for zero eigenvalue is", tolinv, "\n"))
    for (i in 1:length(eigW$values)) {
        if (Ec[i] < tolinv) {
            E[i] <- Ec[i] <- Ec2[i] <- 0
        } else {
            E[i] <- sqrt(1/E[i])
            Ec[i] <- sqrt(1/Ec[i])
            Ec2[i] <- (1/Ec2[i])
        }
    }
    invcW <- diag(Ec)
    irE <- diag(E)
    ZtZ <- irE %*% t(U) %*% t(X) %*% X %*% U %*% irE
    eigZ <- eigen(ZtZ,symmetric=TRUE)
    useEig <- min((ng-1), l)
    A <- Re(eigZ$vectors[, 1:useEig])
    CV <- U %*% invcW %*% A
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
