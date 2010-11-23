mc.unify.list<-function(surpath="sur",lm.data,ref,ext=".ply")
{	
	n<-dim(lm.data)[[3]]
	lmnames<-dimnames(lm.data)[[3]]
	uni.list<-list()
	proclist<-list()
	meshlist<-list()
	proc<-mc.procGPA(lm.data,CSinit=TRUE)
	#nam<-names(meshlist)
	#namref<-dimnames(proc)[[3]][ref]
	#jref<-grep(namref,nam)
	#njref<-names(meshlist)[[jref]]
	
	### read meshes from file ###	
	for (i in 1:n)
		{proclist[[i]]<-file2mesh(paste("sur/",lmnames[i],ext,sep=""))
		
	### create superimposed meshes ### 
		meshlist[[i]]<-rotmesh.onto(proclist[[i]],lm.data[,,i],proc$rotated[,,i],scale=TRUE)$mesh
		}
	
	for (i in 1:n)
	{	#nm0<-dimnames(proc)[[3]][i]
		#j<-grep(nm0,nam)
		#nj<-names(meshlist)[[j]]
		if (i == ref)
			{uni.list[[i]]<-meshlist[[i]]
			
			}
		else
			{uni.list[[i]]<-mc.unify.mesh(meshlist[[ref]],meshlist[[i]],proc$rotated[,,ref],proc$rotated[,,i])$unimesh
			}
		
		
		
		
	}
		names(uni.list)<-lmnames
		names(proclist)<-lmnames
		return(list(uni.list=uni.list,proclist=proclist))
}

