groupPCA <- function(dataarray, groups, rounds = 10000,tol=1e-10,cv=TRUE,mc.cores=detectCores())
  {
    registerDoParallel(cores=mc.cores)### register parallel backend
        
    pmatrix.proc <- NULL
    proc.distout <- NULL
    
    lev<-NULL	
    if (is.character(groups))
      {
        groups<-as.factor(groups)
      }
    if (is.factor(groups))
      {
        factors <- groups
        lev<-levels(groups)
        levn<-length(lev)
        group<-list()
        count<-1
        groupcheck<-0
        for (i in 1:levn)
          {
            tmp0 <- which(groups==lev[i])	
            if (length(tmp0) != 0)
              {			
                group[[count]]<-tmp0
                groupcheck[count]<-i
                count<-count+1
                
              }
            if (length(tmp0)==1)
                {
                  cv=FALSE
                  warning("group with one entry found - crossvalidation will be disabled.")
               }
           
          }
        lev<-lev[groupcheck]
        groups<-group
      }
    
    N <- dataarray
    b <- groups
    
    if (length(dim(N)) == 3) 
		{ n <- dim(N)[3]
                  k <- dim(N)[1]
                  m <- dim(N)[2]
                  l <- k * m
                  ng <- length(groups)
                  if (length(unlist(groups)) != n)
                    {
                      warning("group affinity and sample size not corresponding!")
                    }
                  
                  nwg <- c(rep(0, ng))
                  for (i in 1:ng) 
                    {nwg[i] <- length(b[[i]])
                   }
                  
                  B <- matrix(0, n, m * k)
                  for (i in 1:n) 
                    {
                      B[i, ] <- as.vector(N[, , i])
                    }
                  
                  Gmeans <- matrix(0, ng, m * k)
                  for (i in 1:ng)
                    {
                      Gmeans[i, ] <- as.vector(apply(N[, , b[[i]]], c(1:2),mean))
                    }
                  Grandm <- as.vector(apply(Gmeans, c(1:2), mean))
                  Tmatrix<-B
                  B<-t(t(B)-Grandm)
                  Amatrix <- B
                }
    else
      {
        n <- dim(N)[1]
        l <- dim(N)[2]
        if (length(unlist(groups)) != n)
          {warning("group affinity and sample size not corresponding!")
         }
        ng <- length(groups)
        nwg <- c(rep(0, ng))
        for (i in 1:ng) {
          nwg[i] <- length(b[[i]])
        }
        B <- as.matrix(N)
                                        #Amatrix <- B
        Gmeans <- matrix(0, ng, l)
        for (i in 1:ng)
          {
            Gmeans[i, ] <- apply(N[b[[i]], ], 2, mean)
          }
        Grandm <- apply(Gmeans, 2, mean)
        Tmatrix <- B
	B<-t(t(B)-Grandm)
	Amatrix <- B
      }
    resB <- (Gmeans - (c(rep(1, ng)) %*% t(Grandm)))
    eigenGmeans <- eigen(cov(resB))
    
    valScores <- which(eigenGmeans$values > tol)
    groupScores <- B%*%(eigenGmeans$vectors[,valScores])
    groupPCs=eigenGmeans$vectors[,valScores]
    
    ###### create a neat variance table for the groupmean PCA ###### 
    values <- eigenGmeans$values[valScores]
    if (length(values)==1)
          {Var<-values}
        else
          {
           Var<-matrix(NA,length(values),3)
           Var[,1]<-values
            
            for (i in 1:length(values))
              {
                Var[i,2]<-(values[i]/sum(values))*100
              }
            Var[1,3]<- Var[1,2]
            for (i in 2:length(values))
              {         
                Var[i,3]<-Var[i,2]+ Var[i-1,3]
              }
            colnames(Var)<-c("eigenvalues","% Variance","Cumulative %")
          }

### calculate between group distances ###
  
        proc.disto<-matrix(0, ng, ng)
        if(!is.null(lev))
          {
            rownames(proc.disto)<-lev
            colnames(proc.disto)<-lev
          }	
        for (j1 in 1:(ng - 1)) 
          {
            for (j2 in (j1 + 1):ng) 
              {
                proc.disto[j2, j1] <- sqrt(sum((Gmeans[j1, ]- Gmeans[j2,])^2))
              }
          }
    proc.distout<-as.dist(proc.disto)
    
### Permutation Test for Distances	
    if (rounds != 0) 
      {
        pmatrix.proc <- matrix(NA, ng, ng) ### generate distance matrix Euclidean
        if(!is.null(lev))
          {
            rownames(pmatrix.proc)<-lev
            colnames(pmatrix.proc)<-lev
          }
        rounproc<-function(i)
          {
            b1 <- list()
            dist.mat<-matrix(0,ng,ng)
            shake <- sample(1:n)
            Gmeans1 <- matrix(0, ng, l)
            l1 <- 0
            for (j in 1:ng) 
              {
                b1[[j]] <- c(shake[(l1 + 1):(l1 + (length(b[[j]])))])
                l1 <- l1 + length(b[[j]])
                tmpmat <- Amatrix[b1[[j]],]
                if ( length(b[[j]])==1)
                  tmpmat <- t(as.matrix(tmpmat))
                Gmeans1[j, ] <- apply(tmpmat, 2, mean)
              }
            
            for (j1 in 1:(ng - 1)) 
              {
                for (j2 in (j1 + 1):ng) 
                  {
                    dist.mat[j2, j1] <-sqrt(sum((Gmeans1[j1, ]-Gmeans1[j2, ])^2))
                  }
              }
            
            return(dist.mat)
          }
        
        dist.mat.proc<- array(0, dim = c(ng, ng, rounds))	
        a.list <- foreach(i=1:rounds)%dopar%rounproc(i)
       
        for (i in 1:rounds)
          {
            dist.mat.proc[,,i]<-a.list[[i]]
          }
        for (j1 in 1:(ng - 1)) 
          {
            for (j2 in (j1 + 1):ng) 
              {
                sorti <- sort(dist.mat.proc[j2, j1,])
                
                if (max(sorti) < proc.disto[j2, j1]) 
                  {
                    pmatrix.proc[j2, j1] <- 1/rounds
                  }
                else 
                  {
                    marg <- min(which(sorti >= proc.disto[j2, j1]))
                    pmatrix.proc[j2, j1] <- (rounds - (marg-1))/rounds
                  }
              }
          }
        pmatrix.proc<-as.dist(pmatrix.proc)
      }

    crovafun <- function(x)
      {
        crovtmp <- groupPCAcrova(Tmatrix[-x,],factors[-x],tol=tol,groupPCs=groupPCs)
        out <- (Tmatrix[x,]-crovtmp$Grandmean)%*%crovtmp$PCs      
        return(out)
      }
    CV=NULL
    if (cv)
      {
        crossval <- foreach(i=1:n) %dopar% crovafun(i)
        
        CV <- groupScores
        for (i in 1:n)
          {
            CV[i,] <- crossval[[i]]
          }
      }
     
    return(list(eigenvalues=values,groupPCs=eigenGmeans$vectors[,valScores],Variance=Var,Scores=groupScores,probs=pmatrix.proc,groupdists=proc.distout,groupmeans=Gmeans,Grandmean=Grandm,CV=CV))
        
  }
