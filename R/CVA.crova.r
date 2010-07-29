CVA.crova<-function(dataarray,groups,weighting=TRUE,tolinv=1e-10,ind=0)
{   
    
    
    N<-dataarray
    b<-groups
   
        # starting with a matrix (e.g. PC-Scores)
        n<-dim(N)[1] 
        l<-dim(N)[2]
        ng<-length(b)
	nwg<-c(rep(0,ng))
      for (i in 1:ng)
            {nwg[i]<-length(b[[i]])}
        B<-N
        Amatrix<-B
        Gmeans<-matrix(0,ng,l)
    	
      for (i in 1:ng)
          {#print(N[b[[i]],])
		Gmeans[i,]<-apply(N[b[[i]],],2,mean)}
    
      Grandm<-apply(N[-ind,],2,mean)
      
      
    
      resB<-(Gmeans-(c(rep(1,ng))%*%t(Grandm))) # calculate Between groups Cov Matrix
    if (weighting==TRUE)
    {for (i in 1:ng)
      { resB[i,]<-sqrt(nwg[i])*resB[i,]}
        X<-resB   
    }
    else
      {
      X<-sqrt((n-1)/ng)*resB
    }
      
    #for ( i in 1:ng)
    #{ 
    #  B[b[[i]],]<-B[b[[i]],]-(c(rep(1,length(b[[i]])))%*%t(Gmeans[i,]))
     # }
   
    covW<-0
	for (i in 1:ng) {
        covW <- covW + (cov(B[b[[i]],])*(length(b[[i]])-1))
    }
                     
 #calc within groups covariance
    {covW<-covW+(B[i,]%*%t(B[i,]))
    }
    W<-covW
    covW<-covW/(n-1-ng)
    eigW<-eigen(W)
    eigcoW<-eigen(covW)
    U<-eigW$vectors
    E<-eigW$values
    Ec<-eigcoW$values
    Ec2<-Ec
 	if (min(E) < tolinv)
    {  
        for (i in 1:length(eigW$values))
        {
         if (Ec[i]<tolinv)
          { E[i]<-0
            Ec[i]<-0
            Ec2[i]<-0
          }
        else
          { E[i]<-sqrt(1/E[i])
            Ec[i]<-sqrt(1/Ec[i])
            Ec2[i]<-(1/Ec2[i])
          }
      }
     }
    else 
     {
        for (i in 1:length(eigW$values))
        { E[i]<-sqrt(1/E[i])
            Ec[i]<-sqrt(1/Ec[i])
            Ec2[i]<-(1/Ec2[i])
          }
      }
    irE<-diag(E)
    invcW<-diag(Ec)
    ZtZ<-irE%*%t(U)%*%t(X)%*%X%*%U%*%irE
    eigZ<-eigen(ZtZ)
    A<-eigZ$vectors[,1:(ng-1)]
    CV<-U%*%invcW%*%A
    #CVvis<-covW%*%CV
    #CVscores<-Amatrix%*%CV
    roots<-eigZ$values[1:(ng-1)]
    
    
    
  
      
    ##
    return(list(CV=CV))
}
  
  
  
