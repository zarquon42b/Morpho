CVA <- function (dataarray, groups, weighting = TRUE, tolinv = 1e-10,plot = TRUE, rounds = 0, cv = FALSE, mc.cores=detectCores()) 
{
    if(.Platform$OS.type == "windows")
        mc.cores <- 1

    groups <- factor(groups)
    lev <- levels(groups)
    ng <- length(lev)
    gsizes <- as.vector(tapply(groups, groups, length))
    if (1 %in% gsizes) {
        cv <- FALSE
        warning("group with one entry found - crossvalidation will be disabled.")
    }
    N <- dataarray
    if (length(dim(N)) == 3) 
        N <- vecx(N)
    N <- as.matrix(N)
    n <- dim(N)[1]
    l <- dim(N)[2]
    if (length(groups) != n)
        warning("group affinity and sample size not corresponding!")
    
    n3 <- FALSE
    if (length(dim(N)) == 3) {
        k <- dim(N)[1]
        m <- dim(N)[2]
        N <- vecx(N)
        n3 <- TRUE
    }
    if (is.vector(N) || dim(N)[2] == 1)
        stop("data should contain at least 2 variable dimensions")
    
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng) {
        if (gsizes[i] > 1)
            Gmeans[i, ] <- apply(N[groups==lev[i], ], 2, mean)
        else
            Gmeans[i, ] <- N[groups==lev[i], ]
    }
    Grandm <- as.vector(apply(Gmeans, 2, mean))
    Tmatrix <- N
    N <- sweep(N, 2, Grandm) #center data according to Grandmean
    resGmeans <- sweep(Gmeans, 2, Grandm)
    
    if (weighting) {
        for (i in 1:ng) 
            resGmeans[i, ] <- sqrt(gsizes[i]) * resGmeans[i, ]
        X <- resGmeans
    } else 
        X <- sqrt(n/ng) * resGmeans

    
    
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
    eigZ <- eigen(ZtZ,symmetric=TRUE)
    useEig <- min((ng-1), l)
    A <- Re(eigZ$vectors[, 1:useEig])
    CV <- U %*% invcW %*% A
    CVvis <- covW %*% CV
    CVscores <- N %*% CV

    roots <- eigZ$values[1:useEig]
    if (length(roots) == 1) {
        Var <- matrix(roots, 1, 1)
        colnames(Var) <- "Canonical root"
    } else {
        Var <- matrix(NA, length(roots), 3)
        Var[, 1] <- as.vector(roots)
        for (i in 1:length(roots)) 
            Var[i, 2] <- (roots[i]/sum(roots)) * 100
        Var[1, 3] <- Var[1, 2]
        for (i in 2:length(roots))
            Var[i, 3] <- Var[i, 2] + Var[i - 1, 3]
        colnames(Var) <- c("Canonical roots", "% Variance", "Cumulative %")
    }
    if (plot == TRUE && ng == 2) {
        histGroup(CVscores,groups = groups)
    }
    U2 <- eigcoW$vectors
    winv <- U2 %*% (diag(Ec2)) %*% t(U2)
    disto <- matrix(0, ng, ng)
    rownames(disto) <- colnames(disto) <- lev
         
    #proc.distout <- NULL
    #pmatrix <- NULL
    #pmatrix.proc <- NULL

### calculate Mahalanobis Distance between Means
    for (i in 1:(ng - 1)) {
        for (j in (i + 1):ng) 
            disto[j, i] <- sqrt((Gmeans[i, ] - Gmeans[j,]) %*% winv %*% (Gmeans[i, ] - Gmeans[j, ]))
    }

### calculate Procrustes Distance between Means or Euclidean for 
    proc.disto <- matrix(0, ng, ng)

    if (!is.null(lev))
        colnames(proc.disto) <- rownames(proc.disto) <- lev

    for (i in 1:(ng - 1)) {
        for (j in (i + 1):ng) {
            if (n3) 
                proc.disto[j, i] <- angle.calc(Gmeans[i, ], Gmeans[j,])
            else
                proc.disto[j, i] <- sqrt(sum((Gmeans[i, ]- Gmeans[j,])^2))
        }
    }

### Permutation Test for Distances	
    Dist <- .CVAdists(N, Tmatrix, groups, rounds, proc.dist=proc.disto, maha.dist=disto, mc.cores,n3, winv )
    if (n3) {
        Grandm <- matrix(Grandm, k,m)
        groupmeans <- array(as.vector(t(Gmeans)), dim = c(k,m,ng))
    } else 
        groupmeans <- Gmeans

    CVcv <- NULL
    if (cv == TRUE) {
        CVcv <- CVscores
        ng <- length(groups)
        a.list <- as.list(1:n)
        crovafun <- function(i)
            {
                tmp <- .CVAcrova(Tmatrix[-i, ],groups=groups[-i],test=CV, tolinv = tolinv, weighting=weighting)
                out <- (Tmatrix[i, ]-tmp$Grandmean) %*% tmp$CV
                return(out)
            }
        a.list <- mclapply(a.list, crovafun, mc.cores=mc.cores)
        for (i in 1:n)
            CVcv[i,] <- a.list[[i]]
    }
    return(list(CV = CV, CVscores = CVscores, Grandm = Grandm,
                groupmeans = groupmeans, Var = Var, CVvis = CVvis,
                Dist = Dist,CVcv = CVcv
                ))
}
