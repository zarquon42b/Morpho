CVA <- function (dataarray, groups, weighting = TRUE, tolinv = 1e-10,plot = TRUE, rounds = 0, cv = FALSE, mc.cores=detectCores()) 
{
    if(.Platform$OS.type == "windows")
        {
            mc.cores=1
        }
    lev <- NULL
    k <- 1
    m <- 1
    if (is.character(groups))
        groups <- as.factor(groups)
    
    if (is.factor(groups))
        {
            lev <- levels(groups)
            levn <- length(lev)
            group <- list()
            count <- 1
            groupcheck <- 0
            for (i in 1:levn)
                {
                    tmp0 <- which(groups==lev[i])	
                    if (length(tmp0) != 0)
                        {			
                            group[[count]] <- tmp0
                            groupcheck[count] <- i
                            count <- count+1
                        }
                }
            lev <- lev[groupcheck]
            groups <- group
        }
    
    N <- dataarray
    n3 <- FALSE
    if (length(dim(N)) == 3) 
        {
            k <- dim(N)[1]
            m <- dim(N)[2]
            N <- vecx(N)
            n3 <- TRUE
        }
    if (is.vector(N) || dim(N)[2] == 1)
        stop("data should contain at least 2 variable dimensions")
    
    n <- dim(N)[1]
    l <- dim(N)[2]
    if (length(unlist(groups)) != n)
        {
            warning("group affinity and sample size not corresponding!")
        }
    ng <- length(groups)
    nwg <- c(rep(0, ng))
    for (i in 1:ng)
        {
            nwg[i] <- length(groups[[i]])
        }
                                        #Amatrix <- N
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng)
        {
            if (length(groups[[i]]) > 1)
                Gmeans[i, ] <- apply(N[groups[[i]], ], 2, mean)
            else
                Gmeans[i, ] <- N[groups[[i]], ]
        }
    Grandm <- as.vector(apply(Gmeans, 2, mean))
    Tmatrix <- N
    N <- t(t(N)-Grandm)#center data according to Grandmean
    
    resN <- (Gmeans - (c(rep(1, ng)) %*% t(Grandm)))
    
    if (weighting)
        {
            for (i in 1:ng) 
                resN[i, ] <- sqrt(nwg[i]) * resN[i, ]
            X <- resN
        }
    else 
        X <- sqrt(n/ng) * resN
    
    covW <- 0
    for (i in 1:ng)
        {
            if (!is.vector(N[groups[[i]],]))
                covW <- covW + (cov(N[groups[[i]],])*(length(groups[[i]])-1))
            else
                {
                    covW <- covW+diag(1,l)
                    if (cv)
                        {
                            cv <- FALSE
                            warning("group with one entry found - crossvalidation will be disabled.")
                        }
                }
        }
    
    #W <- covW
    eigW <- eigen(covW)
    covW <- covW/(n - ng)
    eigcoW <- eigen(covW)
    U <- eigW$vectors
    E <- eigW$values
    Ec <- eigcoW$values
    Ec2 <- Ec
    
    if (min(E) < tolinv)
        {
            cat(paste("singular Covariance matrix: General inverse is used. Threshold for zero eigenvalue is", 
                      tolinv, "\n"))
            for (i in 1:length(eigW$values)) {
                if (Ec[i] < tolinv) {
                    E[i] <- 0
                    Ec[i] <- 0
                    Ec2[i] <- 0
                }
                else {
                    E[i] <- sqrt(1/E[i])
                    Ec[i] <- sqrt(1/Ec[i])
                    Ec2[i] <- (1/Ec2[i])
                }
            }
        }
    else
        {
            for (i in 1:length(eigW$values))
                {
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
    if (length(roots) == 1)
        {
            Var <- matrix(roots, 1, 1)
            colnames(Var) <- "Canonical root"
        }
    else
        {
            Var <- matrix(NA, length(roots), 3)
            Var[, 1] <- as.vector(roots)
            for (i in 1:length(roots)) 
                Var[i, 2] <- (roots[i]/sum(roots)) * 100
                        Var[1, 3] <- Var[1, 2]
            for (i in 2:length(roots))
                Var[i, 3] <- Var[i, 2] + Var[i - 1, 3]
            colnames(Var) <- c("Canonical roots", "% Variance", "Cumulative %")
        }
    if (plot == TRUE && ng == 2)
        {
            lim <- range(CVscores[,1])+c(-1,1)
            yli <- c(0,0.7)
            coli <- rainbow(2, alpha = 0.5)
            histo <- hist(CVscores,plot=F)
            hist(CVscores[groups[[1]], ], col = coli[1],add=F,breaks=histo$breaks, main = "CVA", xlab = "CV Scores")
            hist(CVscores[groups[[2]], ], col = coli[2], add = TRUE,breaks=histo$breaks)
        }
    U2 <- eigcoW$vectors
    winv <- U2 %*% (diag(Ec2)) %*% t(U2)
    disto <- matrix(0, ng, ng)
    if (!is.null(lev))
        {
            rownames(disto) <- lev
            colnames(disto) <- lev
        }
    proc.distout <- NULL
    pmatrix <- NULL
    pmatrix.proc <- NULL
    
### calculate Mahalanobis Distance between Means
    for (j1 in 1:(ng - 1)) 
        {
            for (j2 in (j1 + 1):ng) 
                disto[j2, j1] <- sqrt((Gmeans[j1, ] - Gmeans[j2,]) %*% winv %*% (Gmeans[j1, ] - Gmeans[j2, ]))
        }
    
### calculate Procrustes Distance between Means or Euclidean for 
    proc.disto <- matrix(0, ng, ng)

    if (!is.null(lev))
        {
            rownames(proc.disto) <- lev
            colnames(proc.disto) <- lev
        }	
    for (j1 in 1:(ng - 1)) 
        for (j2 in (j1 + 1):ng)
            {
                if (n3) 
                    proc.disto[j2, j1] <- angle.calc(Gmeans[j1, ], Gmeans[j2,])
                else
                    proc.disto[j2, j1] <- sqrt(sum((Gmeans[j1, ]- Gmeans[j2,])^2))
            }
    
    
### Permutation Test for Distances	
   # if (rounds != 0) 
        Dist <- .CVAdists(N=N,Tmatrix=Tmatrix, lev=lev, rounds=rounds, proc.dist=proc.disto,maha.dist=disto, mc.cores=mc.cores,n3=n3, winv=winv,groups=groups)
    
    if (n3) 
        {
            Grandm <- matrix(Grandm, k,m)
            groupmeans <- array(as.vector(t(Gmeans)), dim = c(k,m,ng))
        }
    else 
        groupmeans <- Gmeans
    
    CVcv <- NULL
    if (cv == TRUE) {
        CVcv <- CVscores
        ng <- length(groups)
        a.list <- as.list(1:n)
        crova <- function(i3)
            {
                bb <- groups
                for (j in 1:ng) 
                    if (i3 %in% bb[[j]] == TRUE) 
                        bb[[j]] <- bb[[j]][-(which(bb[[j]] == i3))]
                tmp <- .CVAcrova(N,test=CV, bb, tolinv = tolinv,ind=i3,weighting=weighting)
                out <- (N[i3, ]-tmp$Grandmean) %*% tmp$CV
                return(out)
            }
        a.list <- mclapply(a.list,crova,mc.cores=mc.cores)
        for (i in 1:n)
            CVcv[i,] <- a.list[[i]]
    }
    return(list(CV = CV, CVscores = CVscores, Grandm = Grandm,
                groupmeans = groupmeans, Var = Var, CVvis = CVvis,
                Dist = Dist,CVcv = CVcv
                ))
}
