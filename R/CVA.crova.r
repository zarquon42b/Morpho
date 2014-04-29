.CVAcrova <- function(dataarray,groups,test ,weighting=TRUE,tolinv=1e-10)
{   
    groups <- factor(groups)
    N <- dataarray
    n <- dim(N)[1]
    l <- dim(N)[2]
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))                                 
   
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng) {
        if (gsizes[i] > 1)
            Gmeans[i, ] <- apply(N[groups==lev[i], ], 2, mean)
        else
            Gmeans[i, ] <- N[groups==lev[i], ]
    }
    if (weighting) {
        Grandm <- apply(Gmeans*gsizes,2,sum)/n ## calculate weighted Grandmean (thanks to Anne-Beatrice Dufour for the bug-fix)
    } else {
        Grandm <- as.vector(apply(Gmeans, 2, mean))
    }
    N <- sweep(N, 2, Grandm)
    resGmeans <- sweep(Gmeans, 2, Grandm)
    if (weighting) {
        for (i in 1:ng) 
            resGmeans[i, ] <- sqrt(gsizes[i]) * resGmeans[i, ]
        X <- resGmeans
    } else {
        X <- sqrt(n/ng) * resGmeans
    }
    covW <- covW(N, groups)
    eigW <- eigen(covW*(n - ng))
    eigcoW <- eigen(covW)
    U <- eigW$vectors
    E <- eigW$values
    Ec <- eigcoW$values
    Ec2 <- Ec

    if (min(E) < tolinv)
        cat(paste("singular Covariance matrix: General inverse is used. Threshold for zero eigenvalue is", tolinv, "\n"))
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
    eigZ <- svd(ZtZ)
    useEig <- min((ng-1), l)     
    A <- Re(eigZ$v[, 1:useEig])
    CV <- U %*% invcW %*% A
    di <- dim(CV)[2]
    
    for (i in 1:di) {
        rho <- angle.calc(test[,i],CV[,i])
        if (rho > pi/2)
            CV[,i] <- -CV[,i]		
    }
    return(list(CV=CV,Grandmean=Grandm))
}



