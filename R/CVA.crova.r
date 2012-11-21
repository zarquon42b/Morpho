CVA.crova<-function(dataarray,groups,test=test,weighting=TRUE,tolinv=1e-10,ind=0)
{   
    
    	N <- dataarray
    	b <- groups
    	n <- dim(N)[1]
        l <- dim(N)[2]
	#if (length(unlist(groups)) != (n-1))
	#	{warning("group affinity and sample size not corresponding!")
	#	}
        ng <- length(groups)
        nwg <- c(rep(0, ng))
        for (i in 1:ng) {
            nwg[i] <- length(b[[i]])
        }
        B <- as.matrix(N)
        Amatrix <- B
        Gmeans <- matrix(0, ng, l)
         for (i in 1:ng)
           {
             if(nwg[i] > 1)
               Gmeans[i, ] <- apply(N[b[[i]], ], 2, mean)
             else
               Gmeans[i, ] <- N[b[[i]], ]
           }  
        Grandm <- apply(Gmeans, 2, mean)
    
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
       if (!is.vector(B[b[[i]],]))
        covW <- covW + (cov(B[b[[i]],])*(length(b[[i]])-1))
       else
         covW <- covW+diag(1,l)
              }
    W <- covW
    covW <- covW/(n - ng)
    eigW <- eigen(W)
    eigcoW <- eigen(covW)
    U <- eigW$vectors
    E <- eigW$values
    Ec <- eigcoW$values
    Ec2 <- Ec

   if (min(E) < tolinv) {
        #cat(paste("singular Covariance matrix: General inverse is used. Threshold for zero eigenvalue is",             tolinv, "\n"))
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
        #print(ZtZ)
        tt <- try(eigZ <-svd(ZtZ))
        
    A <- Re(eigZ$v[, 1:(ng - 1)])
	CV <- U %*% invcW %*% A
    di<-dim(CV)[2]
	
	for (i in 1:di)
          {
            rho<-angle.calc(test[,i ],CV[,i])
            if (rho$rho > pi/2)
              
              {
                CV[,i]<-CV[,i]*(-1)			
              }
            
          }
        
    return(list(CV=CV,Grandmean=Grandm))
}
  
  
  
