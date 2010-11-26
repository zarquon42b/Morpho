mc.unify.list<-function(lm.data,ref,surpath="sur",ext=".ply",files=NULL,ray=FALSE,tol=NULL)
{	
	n<-dim(lm.data)[[3]]
	lmnames<-dimnames(lm.data)[[3]]
	uni.list<-list()
	proclist<-list()
	meshlist<-list()
	proc<-mc.procGPA(lm.data,CSinit=TRUE)
	warplist<-list()
	
	### read meshes from file ###
	for (i in 1:n)
		{
		if (is.null(files))
			{meshlist[[i]]<-file2mesh(paste("sur/",lmnames[i],ext,sep=""))
			}
		else	{meshlist[[i]]<-file2mesh(paste("sur/",files[i],sep=""))
			}
		
	### create superimposed meshes ### 
		proclist[[i]]<-rotmesh.onto(meshlist[[i]],lm.data[,,i],proc$rotated[,,i],scale=TRUE)$mesh
		}
	
	for (i in 1:n)
	{	#nm0<-dimnames(proc)[[3]][i]
		#j<-grep(nm0,nam)
		#nj<-names(meshlist)[[j]]
		if (i == ref)
			{uni.list[[i]]<-proclist[[i]]
			warplist[[i]]<-proclist[[i]]
			
			}
		else
			
			{tmp<-mc.unify.mesh(proclist[[ref]],proclist[[i]],proc$rotated[,,ref],proc$rotated[,,i],ray=ray,tol=tol)
			uni.list[[i]]<-tmp$unimesh
			warplist[[i]]<-tmp$warp.mesh
			}
		
		
		
		
	}
		names(uni.list)<-lmnames
		names(proclist)<-lmnames
		return(list(uni.list=uni.list,proclist=proclist,warplist=warplist))
}

