mc.permuvec<-function(data,groups,subgroups,rounds=10000,distabs=FALSE)

{	### define groups ####
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
		for (i in 1:levn)
			{	tmp0<-which(groups==lev[i])	
					if (length(tmp0) != 0)
					{			
					group[[count]]<-tmp0
					count<-count+1
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
        for (i in 1:ng) {
            nwg[i] <- length(b[[i]])
        }
        B <- as.matrix(N)
                
    }
	Gmeans <- matrix(0, ng, l) ### calculate mean of subgroup means for all groups
        	for (i in 1:ng) 
				{for (j in 1:nsub)	
					{tmp<-subgroups[[j]][which(subgroups[[j]] %in% groups[[i]])]
					Gmeans[i,]<-Gmeans[i,]+apply(B[tmp,],2,mean)
					}
					Gmeans[i,]<-Gmeans[i,]/nsub
				}
	
	### create empty subgroup mean matrices 
	meanvec<-matrix(NA,ng,l)
	for (i in 1:ng)
		{meanlist[[i]]<-matrix(NA,nsub,l)
		}
	
	### correct for groupmeans ###
	for (i in 1:ng)
		{gn<-length(groups[[i]])
		delt<-matrix(Gmeans[i,],gn,l,byrow=T)
		B[groups[[i]],]<-B[groups[[i]],]-delt
		}

	### calc subgroup means and residual vectors ###
	for (i in 1:ng)
		{for (j in 1:nsub)	
			{tmp<-subgroups[[j]][which(subgroups[[j]] %in% groups[[i]])]
			meanlist[[i]][j,]<-apply(B[tmp,],2,mean)
			}
		meanvec[i,]<-meanlist[[i]][1,]-meanlist[[i]][2,]
		}
	### calc angle compare vector lengths ###
	if (distabs)
		{disto<-abs(sqrt(sum(meanvec[1,]^2))-sqrt(sum(meanvec[2,]^2)))	
		}
	else
		{disto<-sqrt(sum(meanvec[1,]^2))-sqrt(sum(meanvec[2,]^2))
		}
	
	out<-angle.calc(meanvec[1,],meanvec[2,])$rho
	
	### permutate over groups ###
	

	alist<-as.list(1:rounds)	
	testvec<-0
	permuta<-function(x)
	##tt<-foreach(i = 1:rounds) %dopar%
		{tmplist<-list()
		for (i in 1:ng)
			{tmplist[[i]]<-matrix(NA,nsub,l)
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
		for (i in 1:ng)
			{for (j in 1:nsub)	
				{tmp<-subgroups[[j]][which(subgroups[[j]] %in% b1[[i]])]
				tmplist[[i]][j,]<-apply(Btmp[tmp,],2,mean)
				}
				meanvectmp[i,]<-tmplist[[i]][1,]-tmplist[[i]][2,]
			}
		if (distabs)	
			{dist<-dist<-abs(sqrt(sum(meanvectmp[1,]^2))-sqrt(sum(meanvectmp[2,]^2)))
			}	
		else
			{dist<-dist<-sqrt(sum(meanvectmp[1,]^2))-sqrt(sum(meanvectmp[2,]^2))
			}
		return(c(angle.calc(meanvectmp[1,],meanvectmp[2,])$rho,dist))
		}
		
		tt<-mclapply(alist,permuta)
		uns<-unlist(tt)
		angs<-(1:rounds)*2-1
		dists<-uns[2*(1:rounds)]
		sortdist<-sort(dists)
			if (max(sortdist) < disto)
			{probadist <-1/rounds
			}
		else
			{marg<-min(which(sortdist >= disto))
			probadist<-(rounds-marg)/rounds
			}	
		sortang<-sort(uns[angs])
		if (max(sortang) < out)
			{proba <-1/rounds
			}
		else
			{marg<-min(which(sortang >= out))
			proba<-(rounds-marg)/rounds
			}
		
		
		


	return(list(angle=out,dist=disto,meanvec=meanvec,subgroups=subgroups,permutangles=sortang,permudists=sortdist,pangle=proba,pdist=probadist,uns=uns))
	


	}
	
