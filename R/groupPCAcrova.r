groupPCAcrova <- function(dataarray, groups,tol=1e-10,groupPCs,weighting=weighting)
  { 
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
            if (length(tmp0)==0)
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
      N <- vecx(N)
    
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
    Gmeans <- matrix(0, ng, l)
    for (i in 1:ng)
      {
        if(nwg[i] > 1)
          Gmeans[i, ] <- apply(N[b[[i]], ], 2, mean)
        else
           Gmeans[i, ] <- N[b[[i]], ]
      }  
    if (weighting)
      wt <- nwg
    else
      wt <- rep(1,ng)
    wcov <- cov.wt(Gmeans,wt=wt)
    Grandm <- wcov$center
    eigenGmeans <- eigen(wcov$cov)
    resB <- (Gmeans - (c(rep(1, ng)) %*% t(Grandm)))
    Tmatrix <- B
    B <- t(t(B)-Grandm)
    Amatrix <- B
    valScores <- which(eigenGmeans$values > tol)
    groupScores <- B%*%(eigenGmeans$vectors[,valScores])
    PCs <- as.matrix(eigenGmeans$vectors[,valScores])
    
    if (is.null(dim(groupPCs)))
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
