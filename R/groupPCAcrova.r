groupPCAcrova <- function(dataarray, groups,tol=1e-10,groupPCs,weighting=weighting)
  {
   
    pmatrix.proc <- NULL
    proc.distout <- NULL
    
    lev<-NULL	
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
	B<-t(t(B)-Grandm)
	Amatrix <- B
      }
    resB <- (Gmeans - (c(rep(1, ng)) %*% t(Grandm)))
     if (weighting==TRUE)
      {
        tmpcov <- tcrossprod(Gmeans[1,])*0
        for( i in 1:ng)
          {
            tmpcov <- tmpcov+tcrossprod(resB[i,])
          }
        eigenGmeans <- eigen(tmpcov)
      }
    else
      eigenGmeans <- eigen(cov(resB))
   
    valScores <- which(eigenGmeans$values > tol)
    groupScores <- B%*%(eigenGmeans$vectors[,valScores])
    PCs <-(eigenGmeans$vectors[,valScores])
     if (is.null(dim(PCs)))
         {
           PCs <- matrix(PCs,length(PCs),1)
           groupPCs <-  matrix(groupPCs,length(groupPCs),1)
         }
         di<-dim(PCs)[2]
     
	
	for (i in 1:di)
		{
                  rho<-angle.calc(groupPCs[,i ],PCs[,i])
                  if (rho$rho > pi/2)
                    
			{
                          PCs[,i]<-PCs[,i]*(-1)
                        }

		}
   
   
    return(list(PCs=PCs,Grandmean=Grandm))
        
  }
