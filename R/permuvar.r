permuvar<-function(data,groups,rounds)
{	N <- data
    

	if (!is.factor(groups))
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
		b <- groups
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
		distlist<-list()
        B <- as.matrix(N)
		Gmeans <- matrix(0, ng, l)
        for (i in 1:ng) 
			{
 			Gmeans[i, ] <- as.vector(apply(N[ b[[i]],], 2, mean))
			delt<-matrix(Gmeans[i,],nwg[i],l,byrow=T)
			B[groups[[i]],]<-B[groups[[i]],]-delt
			distlist[[i]]<-0
			for (j in 1:length(groups[[i]]))
				{distlist[[i]]<-distlist[[i]]+(sum(B[groups[[i]][j],]^2))
				}
			distlist[[i]]<-distlist[[i]]/nwg[i]
			
			}
			realdist<-(distlist[[1]]-distlist[[2]])
	alist<-0
	for (tt in 1:rounds)
		{	distlist0<-list()
			Gmeans0 <- matrix(0, ng, l)
        	b0<-list()
			shake<-sample(1:n)
			l1 <- 0
			B1<-B
		for (i in 1:ng) 
			{b0[[i]]<-c(shake[(l1 + 1):(l1 + (length(b[[i]])))])
            l1 <- l1 + length(b[[i]])
			
			Gmeans0[i, ] <- as.vector(apply(N[ b0[[i]],], 2, mean))
			delt<-matrix(Gmeans0[i,],nwg[i],l,byrow=T)
			B1[b0[[i]],]<-B1[b0[[i]],]-delt
			distlist0[[i]]<-0
			for (j in 1:length(groups[[i]]))
				{distlist0[[i]]<-distlist0[[i]]+(sum(B1[b0[[i]][j],]^2))
				}
			distlist0[[i]]<-distlist0[[i]]/nwg[i]

			}
		alist[tt]<-abs(distlist0[[1]]-distlist0[[2]])
		}
		
	return(list(alist,realdist))
	


}

