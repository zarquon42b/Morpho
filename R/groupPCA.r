groupPCA <- function(dataarray, groups, rounds = 10000,tol=1e-10,cv=TRUE,mc.cores=detectCores(), weighting=TRUE)
{
    win <- FALSE
    if(.Platform$OS.type == "windows")
        win <- TRUE
    else
        registerDoParallel(cores=mc.cores)### register parallel backend
    
    pmatrix.proc <- NULL
    proc.distout <- NULL
    lev <- NULL	
    
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
    
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng) {
        if(gsizes[i] > 1)
            Gmeans[i, ] <- apply(N[groups==lev[i], ], 2, mean)
        else
            Gmeans[i, ] <- N[groups==lev[i], ]
    }
    if (weighting==TRUE)
        wt <- gsizes
    else
        wt <- rep(1,ng)
    wcov <- cov.wt(Gmeans,wt=wt)
    Grandm <- wcov$center
    eigenGmeans <- eigen(wcov$cov)
                                        #resGmeans <- sweep(Gmeans, 2, Grandm)
    Tmatrix <- N
    N <- sweep(N, 2, Grandm)
    valScores <- which(eigenGmeans$values > tol)
    groupScores <- N%*%(eigenGmeans$vectors[,valScores])
    groupPCs <- eigenGmeans$vectors[,valScores]

###### create a neat variance table for the groupmean PCA ###### 
    values <- eigenGmeans$values[valScores]
    if (length(values) == 1) {
        Var <- values
    } else {
        Var <- matrix(NA,length(values),3)
        Var[,1] <- values
        for (i in 1:length(values)) 
            Var[i,2] <- (values[i]/sum(values))*100
        Var[1,3] <- Var[1,2]
        for (i in 2:length(values))
            Var[i,3] <- Var[i,2]+ Var[i-1,3]
        colnames(Var) <- c("eigenvalues","% Variance","Cumulative %")
    }
### calculate between group distances ###
    proc.disto <- matrix(0, ng, ng)
    if(!is.null(lev)) {
        rownames(proc.disto) <- lev
        colnames(proc.disto) <- lev
    }	
    for (j1 in 1:(ng - 1)) 
        for (j2 in (j1 + 1):ng) 
            proc.disto[j2, j1] <- sqrt(sum((Gmeans[j1, ]- Gmeans[j2,])^2))
    
    proc.distout <- as.dist(proc.disto)

### Permutation Test for Distances	
    if (rounds > 0) {
        pmatrix.proc <- matrix(NA, ng, ng) ### generate distance matrix Euclidean
        if(!is.null(lev)) {
            rownames(pmatrix.proc) <- lev
            colnames(pmatrix.proc) <- lev
        }
        rounproc <- function(i)
            {
                shake <- sample(groups)
                dist.mat <- matrix(0,ng,ng)
                Gmeans.tmp <- matrix(0, ng, l)
                for (j in 1:ng) {
                    if(gsizes[j] > 1)
                        Gmeans.tmp[j, ] <- apply(N[shake==lev[j], ], 2, mean)
                    else
                        Gmeans.tmp[j, ] <- N[shake==lev[j], ]
                }
                dist.mat <- as.matrix(dist(Gmeans.tmp))
                return(dist.mat)
            }
        
        dist.mat.proc <- array(0, dim = c(ng, ng, rounds))
        if(win)
            a.list <- foreach(i=1:rounds)%do%rounproc(i)
        else
            a.list <- foreach(i=1:rounds)%dopar%rounproc(i)
        
        for (i in 1:rounds)
            dist.mat.proc[,,i] <- a.list[[i]]
        
        for (j1 in 1:(ng - 1)) {
            for (j2 in (j1 + 1):ng) {
                sorti <- sort(dist.mat.proc[j2, j1,])
                if (max(sorti) < proc.disto[j2, j1]) {
                    pmatrix.proc[j2, j1] <- 1/rounds
                } else {
                    marg <- min(which(sorti >= proc.disto[j2, j1]))
                    pmatrix.proc[j2, j1] <- (rounds - (marg-1))/rounds
                }
            }
        }
        pmatrix.proc <- as.dist(pmatrix.proc)
    }
    crovafun <- function(x)
        {
            crovtmp <- .groupPCAcrova(Tmatrix[-x,],groups[-x],tol=tol,groupPCs=groupPCs,weighting=weighting)
            out <- as.vector(Tmatrix[x,]-crovtmp$Grandmean) %*% as.matrix(crovtmp$PCs)
            return(out)
        }
    
    CV=NULL
    if (cv) {
        if (win)
            crossval <- foreach(i=1:n) %do% crovafun(i)
        else
            crossval <- foreach(i = 1:n) %dopar% crovafun(i)
        CV <- groupScores
        for (i in 1:n) {
            if (is.matrix(CV))
                CV[i,] <- crossval[[i]]
            else
                CV[i] <- crossval[[i]]
        }
    }
    return(list(eigenvalues=values,groupPCs=eigenGmeans$vectors[,valScores],Variance=Var,Scores=groupScores,probs=pmatrix.proc,groupdists=proc.distout,groupmeans=Gmeans,Grandmean=Grandm,CV=CV))
}
