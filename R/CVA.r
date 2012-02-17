CVA <- mc.CVA<-function (dataarray, groups, weighting = TRUE, tolinv = 1e-10,plot = TRUE, rounds = 10000, cv = TRUE) 
{
  
  lev<-NULL
  if(.Platform$OS.type == "windows")
     {
     # registerDoParallel(cores=cores)
      tmplapply <- lapply
    }
    else
      {
      #  registerDoMC(cores=cores)
        tmplapply <- mclapply
      }	
  if (is.character(groups))
    {
      groups<-as.factor(groups)
    }
  if (is.factor(groups))
    {
      lev<-levels(groups)
      levn<-length(lev)
      group<-list()
      count<-1
      groupcheck<-0
      for (i in 1:levn)
        {
          tmp0<-which(groups==lev[i])	
          if (length(tmp0) != 0)
            {			
              group[[count]]<-tmp0
              groupcheck[count]<-i
              count<-count+1
              
            }
        }
      lev<-lev[groupcheck]
      groups<-group
    }
  
  N <- dataarray
  b <- groups
  
  if (length(dim(N)) == 3) 
    {
      n <- dim(N)[3]
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
        {
          nwg[i] <- length(b[[i]])
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
      Grandm <- as.vector(apply(Gmeans, 2, mean))
      Tmatrix<-B
      B<-t(t(B)-Grandm)
      Amatrix <- B
    }
  else
    {
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
          nwg[i] <- length(b[[i]])
        }
      B <- as.matrix(N)
                                        #Amatrix <- B
      Gmeans <- matrix(0, ng, l)
      for (i in 1:ng) {
        Gmeans[i, ] <- apply(N[b[[i]], ], 2, mean)
      }
       Grandm <- as.vector(apply(Gmeans, 2, mean))
      Tmatrix<-B
      B<-t(t(B)-Grandm)
      Amatrix <- B
    }
  resB <- (Gmeans - (c(rep(1, ng)) %*% t(Grandm)))
  
  if (weighting == TRUE) {
    for (i in 1:ng) {
      resB[i, ] <- sqrt(nwg[i]) * resB[i, ]
    }
    X <- resB
  }
  else {
    X <- sqrt(n/ng) * resB
  }
  
  covW <- 0
  for (i in 1:ng) {
    covW <- covW + (cov(B[b[[i]],])*(length(b[[i]])-1))
  }
  
  
  W <- covW
  covW <- covW/(n - ng)
  eigW <- eigen(W)
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
  else {
    for (i in 1:length(eigW$values)) {
      E[i] <- sqrt(1/E[i])
      Ec[i] <- sqrt(1/Ec[i])
      Ec2[i] <- (1/Ec2[i])
    }
  }
  
  
  invcW <- diag(Ec)
  irE <- diag(E)
  ZtZ <- irE %*% t(U) %*% t(X) %*% X %*% U %*% irE
  eigZ <- eigen(ZtZ,symmetric=TRUE)
  A <- Re(eigZ$vectors[, 1:(ng - 1)])
  CV <- U %*% invcW %*% A
  CVvis <- covW %*% CV
  CVscores <- Amatrix %*% CV
  
  roots <- eigZ$values[1:(ng - 1)]
  if (length(roots) == 1) {
    Var <- matrix(roots, 1, 1)
    colnames(Var) <- "Canonical root"
  }
  else {
    Var <- matrix(NA, length(roots), 3)
    Var[, 1] <- as.vector(roots)
    for (i in 1:length(roots)) {
      Var[i, 2] <- (roots[i]/sum(roots)) * 100
    }
    Var[1, 3] <- Var[1, 2]
    for (i in 2:length(roots)) {
      Var[i, 3] <- Var[i, 2] + Var[i - 1, 3]
    }
    colnames(Var) <- c("Canonical roots", "% Variance", "Cumulative %")
  }
  if (plot == TRUE && ng == 2) {
    lim<-range(CVscores[,1])+c(-1,1)
    yli<-c(0,0.7)
    coli <- rainbow(2, alpha = 0.5)
    histo <- hist(CVscores,plot=F)
    hist(CVscores[b[[1]], ], col = coli[1],add=F,breaks=histo$breaks, main = "CVA", xlab = "CV Scores")
    hist(CVscores[b[[2]], ], col = coli[2], add = TRUE,breaks=histo$breaks)
  }
  U2 <- eigcoW$vectors
  winv <- U2 %*% (diag(Ec2)) %*% t(U2)
  disto <- matrix(0, ng, ng)
  if(!is.null(lev))
    {rownames(disto)<-lev
     colnames(disto)<-lev
   }
  proc.distout<-NULL
  pmatrix<-NULL
  pmatrix.proc<-NULL
  
### calculate Mahalanobis Distance between Means
  for (j1 in 1:(ng - 1)) 
    {
      for (j2 in (j1 + 1):ng) 
        {disto[j2, j1] <- sqrt((Gmeans[j1, ] - Gmeans[j2,]) %*% winv %*% (Gmeans[j1, ] - Gmeans[j2, ]))
       }
    }
  
### calculate Procrustes Distance between Means
  if (length(dim(N)) == 3)
    {
      proc.disto<-matrix(0, ng, ng)
      if(!is.null(lev))
        {rownames(proc.disto)<-lev
         colnames(proc.disto)<-lev
       }	
      for (j1 in 1:(ng - 1)) 
        {
          for (j2 in (j1 + 1):ng) 
            {proc.disto[j2, j1] <- angle.calc(Gmeans[j1, ], Gmeans[j2,])$rho
           }
        }
      proc.distout<-as.dist(proc.disto)
      
    }
        else
          {
            proc.disto<-matrix(0, ng, ng)
            if(!is.null(lev))
              {rownames(proc.disto)<-lev
               colnames(proc.disto)<-lev
             }	
            for (j1 in 1:(ng - 1)) 
              {
       		for (j2 in (j1 + 1):ng) 
                  {proc.disto[j2, j1] <- sqrt(sum((Gmeans[j1, ]- Gmeans[j2,])^2))
                 }
              }
            proc.distout<-as.dist(proc.disto)
            
          }
  
### Permutation Test for Distances	
	if (rounds != 0) 
	{pmatrix <- matrix(NA, ng, ng) ### generate distance matrix Mahal
	if(!is.null(lev))
		{rownames(pmatrix)<-lev
		colnames(pmatrix)<-lev
		}
	
	
	
	
	roun<-function(i)
	{b1 <- list(numeric(0))
		dist.mat<-matrix(0,ng,ng)
            shake <- sample(1:n)
            Gmeans1 <- matrix(0, ng, l)
            l1 <- 0
            for (j in 1:ng) {
                b1[[j]] <- c(shake[(l1 + 1):(l1 + (length(b[[j]])))])
                l1 <- l1 + length(b[[j]])
                Gmeans1[j, ] <- apply(Amatrix[b1[[j]], ], 2, 
                  mean)
            }
            for (j1 in 1:(ng - 1)) {
                for (j2 in (j1 + 1):ng) {
                  dist.mat[j2, j1] <- sqrt((Gmeans1[j1, ] - Gmeans1[j2, ]) %*% winv %*% (Gmeans1[j1,] - Gmeans1[j2, ]))
                }
            }
	return(dist.mat)
	}
	
        dist.mat<- array(0, dim = c(ng, ng, rounds))
        a.list<-as.list(1:rounds)
	a.list<-tmplapply(a.list,roun)
	for (i in 1:rounds)
		{dist.mat[,,i]<-a.list[[i]]
		#print(a.list)[[i]]
		}
        #pmatrix <- matrix(0, ng, ng)
	
	
        for (j1 in 1:(ng - 1)) 
		{for (j2 in (j1 + 1):ng) 
			{sorti <- sort(dist.mat[j2, j1, ])
                	if (max(sorti) < disto[j2, j1]) 
				{pmatrix[j2, j1] <- 1/rounds
                		}
               		else 
                  		{marg <- min(which(sorti >= disto[j2, j1]))
                  		pmatrix[j2, j1] <- (rounds - marg+1)/rounds
                		}
            		}
        	}

    	if (length(dim(N)) == 3)
		{
                  pmatrix.proc <- matrix(NA, ng, ng) ### generate distance matrix ProcDist for Landmark configurations
		if(!is.null(lev))
			{rownames(pmatrix.proc)<-lev
			colnames(pmatrix.proc)<-lev
			}
		roun.proc<-function(i)
			{b1 <- list()
			dist.mat<-matrix(0,ng,ng)
            		shake <- sample(1:n)
            		Gmeans1 <- matrix(0, ng, l)
            		l1 <- 0
            		for (j in 1:ng) 
				{b1[[j]] <- c(shake[(l1 + 1):(l1 + (length(b[[j]])))])
                		l1 <- l1 + length(b[[j]])
                		Gmeans1[j, ] <- apply(Tmatrix[b1[[j]], ], 2, mean)
            			}
			
            		for (j1 in 1:(ng - 1)) 
				{for (j2 in (j1 + 1):ng) 
					{dist.mat[j2, j1] <- angle.calc(Gmeans1[j1, ],Gmeans1[j2, ])$rho
                			}
            			}
			
			return(dist.mat)
			}
		
        	dist.mat.proc<- array(0, dim = c(ng, ng, rounds))
        	a.list<-as.list(1:rounds)
		a.list<-tmplapply(a.list,roun.proc)
		
		for (i in 1:rounds)
			{dist.mat.proc[,,i]<-a.list[[i]]
			}
		for (j1 in 1:(ng - 1)) 
			{for (j2 in (j1 + 1):ng) 
				{sorti <- sort(dist.mat.proc[j2, j1, ])
				if (max(sorti) < proc.disto[j2, j1]) 
					{pmatrix.proc[j2, j1] <- 1/rounds
                			}
               			else 
                  			{marg <- min(which(sorti >= proc.disto[j2, j1]))
                  			pmatrix.proc[j2, j1] <- (rounds - marg+1)/rounds
                			}
            			}
        		}
		
		pmatrix.proc<-as.dist(pmatrix.proc)
		}
         else ## case if input is datamatrix only euclidean distances will be calculated
           	{pmatrix.proc <- matrix(NA, ng, ng) ### generate distance matrix ProcDist for Landmark configurations
		if(!is.null(lev))
			{rownames(pmatrix.proc)<-lev
			colnames(pmatrix.proc)<-lev
			}
		roun.proc<-function(i)
			{b1 <- list()
			dist.mat<-matrix(0,ng,ng)
            		shake <- sample(1:n)
            		Gmeans1 <- matrix(0, ng, l)
            		l1 <- 0
            		for (j in 1:ng) 
				{b1[[j]] <- c(shake[(l1 + 1):(l1 + (length(b[[j]])))])
                		l1 <- l1 + length(b[[j]])
                		Gmeans1[j, ] <- apply(Amatrix[b1[[j]], ], 2, mean)
            			}
			
            		for (j1 in 1:(ng - 1)) 
				{for (j2 in (j1 + 1):ng) 
					{dist.mat[j2, j1] <-sqrt(sum((Gmeans1[j1, ]-Gmeans1[j2, ])^2))
                			}
            			}
			
			return(dist.mat)
			}
		
        	dist.mat.proc<- array(0, dim = c(ng, ng, rounds))
        	a.list<-as.list(1:rounds)
		a.list<-tmplapply(a.list,roun.proc)
		
		for (i in 1:rounds)
			{dist.mat.proc[,,i]<-a.list[[i]]
			}
		for (j1 in 1:(ng - 1)) 
			{for (j2 in (j1 + 1):ng) 
				{sorti <- sort(dist.mat.proc[j2, j1, ])
				if (max(sorti) < proc.disto[j2, j1]) 
					{pmatrix.proc[j2, j1] <- 1/rounds
                			}
               			else 
                  			{marg <- min(which(sorti >= proc.disto[j2, j1]))
                  			pmatrix.proc[j2, j1] <- (rounds - marg +1)/rounds
                			}
            			}
        		}
		
		pmatrix.proc<-as.dist(pmatrix.proc)
		}
         
	pmatrix <- as.dist(pmatrix)
	
	}
	
	disto <- as.dist(disto)
	
    	Dist <- list(GroupdistMaha = disto,GroupdistProc=proc.distout, probsMaha = pmatrix,probsProc = pmatrix.proc)
    	if (length(dim(N)) == 3) 
		{Grandm <- matrix(Grandm, k, m)
        	groupmeans <- array(as.vector(t(Gmeans)), dim = c(k,m,ng))
    		}
   
	else 
		{groupmeans <- Gmeans
    		}
    CVcv <- NULL
    if (cv == TRUE) {
        CVcv <- CVscores
        ng <- length(groups)
	a.list<-as.list(1:n)
	crova<-function(i3)
		{bb <- groups
            	for (j in 1:ng) 
			{if (i3 %in% bb[[j]] == TRUE) 
				{bb[[j]] <- bb[[j]][-(which(bb[[j]] == i3))]
                		}
            		}
		
            	tmp <- CVA.crova(Amatrix,test=CV, bb, tolinv = tolinv,ind=i3)
            	out <- (Tmatrix[i3, ]-tmp$Grandmean) %*% tmp$CV
		return(out)
        	}
	a.list<-tmplapply(a.list,crova)
	for (i in 1:n)
	CVcv[i,]<-a.list[[i]]
	
    }
    return(list(CV = CV, CVscores = CVscores, Grandm = Grandm, 
        groupmeans = groupmeans, Var = Var, CVvis = CVvis, Dist = Dist,CVcv = CVcv))
}
