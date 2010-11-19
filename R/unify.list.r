unify.list<-function(meshlist,proc,ref)
{	
	n<-dim(proc)[[3]]
	uni.list<-list()
	nam<-names(meshlist)
	namref<-dimnames(proc)[[3]][ref]
	jref<-grep(namref,nam)
	njref<-names(meshlist)[[jref]]
	
	
	for (i in 1:n)
	{	nm0<-dimnames(proc)[[3]][i]
		j<-grep(nm0,nam)
		nj<-names(meshlist)[[j]]
		if (i == ref)
			{uni.list[[i]]<-meshlist[[j]]
			
			}
		else
			{uni.list[[i]]<-unify.mesh(meshlist[[jref]],meshlist[[j]],proc[,,ref],proc[,,i])$unimesh
			}
		
		
		names(uni.list)[i]<-nj
		#print(j)
	}
		
		return(uni.list)
	}

