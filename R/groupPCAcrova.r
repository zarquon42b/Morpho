.groupPCAcrova <- function(dataarray, groups,tol=1e-10,groupPCs,weighting=weighting)
    { 
        lev <- NULL	
        if (is.character(groups))
            {
                groups <- as.factor(groups)
            }
        if (is.factor(groups))
            {
                factors <- groups
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
        if (length(dim(N)) == 3) 
            N <- vecx(N)
        
        n <- dim(N)[1]
        l <- dim(N)[2]
        if (length(unlist(groups)) != n)
            {warning("group affinity and sample size not corresponding!")
         }
        ng <- length(groups)
        nwg <- c(rep(0, ng))
        for (i in 1:ng) {
            nwg[i] <- length(groups[[i]])
        }
        N <- as.matrix(N)
        Gmeans <- matrix(0, ng, l)
        for (i in 1:ng)
            {
                if(nwg[i] > 1)
                    Gmeans[i, ] <- apply(N[groups[[i]], ], 2, mean)
                else
                    Gmeans[i, ] <- N[groups[[i]], ]
            }  
        if (weighting)
            wt <- nwg
        else
            wt <- rep(1,ng)
        wcov <- cov.wt(Gmeans,wt=wt)
        Grandm <- wcov$center
        eigenGmeans <- eigen(wcov$cov)
        resN <- (Gmeans - (c(rep(1, ng)) %*% t(Grandm)))
        Tmatrix <- N
        N <- t(t(N)-Grandm)
        valScores <- which(eigenGmeans$values > tol)
        groupScores <- N%*%(eigenGmeans$vectors[,valScores])
        PCs <- as.matrix(eigenGmeans$vectors[,valScores])
        
        if (is.null(dim(groupPCs)))
            {
                PCs <- matrix(PCs,length(PCs),1)
                groupPCs <- matrix(groupPCs,length(groupPCs),1)
            }
        di <- dim(PCs)[2]
        for (i in 1:di)
            {
                rho <- angle.calc(groupPCs[,i ],PCs[,i])
                if (rho > pi/2)
                    {
                        PCs[,i] <- PCs[,i]*(-1)
                    }
            }
        return(list(PCs=PCs,Grandmean=Grandm))
    }
