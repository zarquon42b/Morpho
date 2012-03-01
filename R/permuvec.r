permuvec <- mc.permuvec<-function(data,groups,subgroups,rounds=10000,scale=TRUE,tol=1e-10)

{
  registerDoParallel() ##register parallel backend
  
### define groups ####
  rawgroup<-groups	
  lev<-NULL	
  if (is.character(groups) || is.numeric(groups))
    {groups<-as.factor(groups)
   }
  if (is.factor(groups))
    {
      lev<-levels(groups)
      levn<-length(lev)
      group<-list()
      count<-1
      groupcheck<-0
      for (i in 1:levn)
        {	tmp0<-which(groups==lev[i])	
                if (length(tmp0) != 0)
                  {			
                    group[[count]]<-tmp0
                    count<-count+1
                    groupcheck[count]<-i
                  }
              }
      groups<-group
    }
  
  levsub<-NULL	
  if (is.character(subgroups) || is.numeric(subgroups))
    {subgroups<-as.factor(subgroups)
   }
  if (is.factor(subgroups))
    {
      levsub<-levels(subgroups)
      levnsub<-length(levsub)
      subgroup<-list()
      count<-1
      for (i in 1:levnsub)
        {	tmp0<-which(subgroups==levsub[i])
                if (length(tmp0) != 0)
                  {subgroup[[i]]<-tmp0
                   count<-count+1
                 }
              }
      subgroups<-subgroup
    }
  b<-groups
  N <- data
  ng <- length(groups)
  nsub<-length(subgroups)
  meanlist<-list()
  
### prepare data if data is an array ###
  
  if (length(dim(N)) == 3) 
    { n <- dim(N)[3]
      k <- dim(N)[1]
      m <- dim(N)[2]
      l <- k * m
                  
      if (length(unlist(groups)) != n)
        {warning("group affinity and sample size not corresponding!")
       }
      
      nwg <- c(rep(0, ng))
                  for (i in 1:ng) 
                    {nwg[i] <- length(b[[i]])
                   }
      
      B <- matrix(0, n, m * k)
      for (i in 1:n) 
        {B[i, ] <- as.vector(N[, , i])
       }  
    }
  else {
    n <- dim(N)[1]
    l <- dim(N)[2]
    if (length(unlist(groups)) != n)
      {warning("group affinity and sample size not corresponding!")
     }
    ng <- length(groups)
    nwg <- c(rep(0, ng))
    for (i in 1:ng)
      {
        nwg[i] <- length(b[[i]])
      }
    B <- as.matrix(N)
    
  }
  Braw<-B
  nws <- c(rep(0, nsub))
  for (i in 1:nsub)
    {
      nws[i] <- length(subgroups[[i]])
    }
  Gmeans <- matrix(0, ng, l) ### calculate mean of subgroup means for all groups ###
  for (i in 1:ng) 
    {
      for (j in 1:nsub)	
        {
          tmp<-subgroups[[j]][which(subgroups[[j]] %in% groups[[i]])]
          Gmeans[i,]<-Gmeans[i,]+apply(B[tmp,],2,mean)
        }
      Gmeans[i,]<-Gmeans[i,]/nsub
    }
  
### create empty subgroup mean matrices 
  meanvec<-matrix(NA,ng,l)
  if (is.factor (rawgroup ) || is.character(rawgroup))
    {
      rownames(meanvec)<-levels(rawgroup)[groupcheck]
    }
  
  for (i in 1:ng)
    {
      meanlist[[i]]<-matrix(NA,nsub,l)
    }
  
### correct for groupmeans ###
  for (i in 1:ng)
    {
      gn<-length(groups[[i]])
      delt<-matrix(Gmeans[i,],gn,l,byrow=T)
      B[groups[[i]],]<-B[groups[[i]],]-delt
    }
  
### calc subgroup means, residual vectors and pooled within group variance ###
	covW <- 0	
  wgroupvar<-NULL
  for (i in 1:ng)
    {	
      for (j in 1:nsub)	
              {
                tmp<-subgroups[[j]][which(subgroups[[j]] %in% groups[[i]])]
                meanlist[[i]][j,]<-apply(B[tmp,],2,mean)
### calc within subgroups Sum of Squares
                if (scale)
                  {
                    covW<-covW+cov(apply(B[tmp,],2,scale,scale=F))*(length(tmp)-1)
                  }
              }
      
### calc pooled groupspecific within subgroups covariance matrix and overall variance ###
      #print(scale)
      
      meanvec[i,]<-(meanlist[[i]][1,]-meanlist[[i]][2,])      
    }
  
  covW<-covW/(n-(ng*nsub))
  if (!scale)
    {
      covW <-  diag(rep(1,dim(B)[2]))
    }
  mahadist<-NULL
### invert covariance matrix
  coinv<-mpinv(covW,tol=tol)
                                        # print(dim(coinv))
  for (i in 1:ng) ## calc Mahalanobisdistance ### 
    {mahadist[i]<-sqrt(meanvec[i,]%*%coinv%*%meanvec[i,])
                                        # print(dim(sqrt(meanvec[i,]%*%coinv%*%meanvec[i,])))
   }
  
### calc angle compare vector lengths ###
  
  disto<-abs(mahadist[1]-mahadist[2])
  out<-Morpho::angle.calc(meanvec[1,],meanvec[2,])$rho
  
  
### permutate over groups ###	
  
  alist<-as.list(1:rounds)	
  testvec<-0
  permuta<-function(x)
    {
      tmplist<-list()
      for (i in 1:ng)
        {
          tmplist[[i]]<-matrix(NA,nsub,l)
        }
      meanvectmp<-matrix(NA,ng,l)
      Btmp<-B
      b1 <- list(numeric(0))
      shake<-sample(1:n)
      Gmeans1 <- matrix(0, ng, l)
      l1 <- 0
      for (j in 1:ng) 
        {
          b1[[j]] <- c(shake[(l1 + 1):(l1 + (length(b[[j]])))])
          l1 <- l1 + length(b[[j]])
        }
                                        #covW0 <- 0
      for (i in 1:ng)
        {	
          for (j in 1:nsub)	
            {
              tmp<-subgroups[[j]][which(subgroups[[j]] %in% b1[[i]])]
              tmplist[[i]][j,]<-apply(Btmp[tmp,],2,mean)
                                        #covW0<-covW0+cov(apply(B[tmp,],2,scale,scale=F))*(length(tmp)-1)
            }
          
          meanvectmp[i,]<-(tmplist[[i]][1,]-tmplist[[i]][2,])
                                        #meanvectmp[i,]
        }
      
      mahadist0<-NULL
      
      for (i in 1:ng) ## calc Mahalanobisdistance ### 
        {
          mahadist0[i]<-sqrt(meanvectmp[i,]%*%coinv%*%meanvectmp[i,])        
        }
                                        #covW0<-covW0/(n-(ng*nsub))
      dist<-abs(mahadist0[1]-mahadist0[2])
      
      return(c(Morpho::angle.calc(meanvectmp[1,],meanvectmp[2,])$rho,dist))
    }
  
  tt <- foreach(i= 1:rounds,.combine=c) %dopar% permuta(i)
    
   # mclapply(alist,permuta)
  uns <- unlist(tt)
  angs <- (1:rounds)*2-1
  dists <- uns[2*(1:rounds)]
  sortdist <- sort(dists)
  
### calc probabilities ####
  if (max(sortdist) < disto)
    {
      probadist <-1/rounds
    }
  else
    {
      marg<-min(which(sortdist >= disto))
      probadist<-(rounds-marg+1)/rounds
    }      
  sortang<-sort(uns[angs])
  if (max(sortang) < out)
    {proba <-1/rounds
   }
  else
    {
      marg<-min(which(sortang >= out))
      proba<-(rounds-marg+1)/rounds
    }      

  return(list(angle=out,dist=disto,meanvec=meanvec,permutangles=sortang,permudists=sortdist,p.angle=proba,p.dist=probadist,subdist=mahadist))
}
	
