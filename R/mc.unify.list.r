mc.unify.list<-function(lm.data,ref,surpath="sur/",ext=".ply",files=NULL,ray=FALSE,tol=NULL,rawout=NULL,uniout=NULL,keep=TRUE)
{	
	if (is.character(rawout))
		
		{	le<-nchar(rawout)
			if (substr(rawout,le,le)!="/")
				{rawout<-paste(rawout,"/",sep="")
				}
			
			a<-dir(rawout)
			if (length(a)==0)
			{dir.create(rawout)
			}
		} 	
	if (is.character(uniout))
		
		{	le<-nchar(uniout)
			if (substr(uniout,le,le)!="/")
				{uniout<-paste(uniout,"/",sep="")
				}
			
			a<-dir(uniout)
			if (length(a)==0)
			{dir.create(uniout)
			}
		} 	
	n<-dim(lm.data)[[3]]
	lmnames<-dimnames(lm.data)[[3]]
	uni.list<-list()
	#proclist<-list()
	meshlist<-list()
	proc<-mc.procGPA(lm.data,CSinit=TRUE)
	#warplist<-list()
	
	### read reference mesh ###

	refmesh<-file2mesh(paste(surpath,lmnames[ref],ext,sep=""))
	refprocmesh<-rotmesh.onto(refmesh,lm.data[,,ref],proc$rotated[,,ref],scale=TRUE)$mesh
	### read meshes from file ###
	for (i in 1:n)
		{
			if (is.null(files))
				{tmpmesh<-file2mesh(paste(surpath,lmnames[i],ext,sep=""))
				}
			else	
				{tmpmesh<-file2mesh(paste(surpath,files[i],sep=""))
				}
		
	### create superimposed meshes ### 
		
		procmesh<-rotmesh.onto(tmpmesh,lm.data[,,i],proc$rotated[,,i],scale=TRUE)$mesh
		
		if (i == ref)
			{
			if (keep)
			{uni.list[[i]]<-refprocmesh
			#warplist[[i]]<-refprocmesh
			#proclist[[i]]<-refprocmesh
			}
			if (is.character(rawout))
				{
				mesh2ply(refmesh,paste(rawout,lmnames[i],sep=""))
				}

			if (is.character(uniout))
				{mesh2ply(refprocmesh,paste(uniout,lmnames[i],sep=""))
				}
			
			}
		else
			
			{tmp<-unify.mesh(refprocmesh,procmesh,proc$rotated[,,ref],proc$rotated[,,i],ray=ray,tol=tol)
			if (keep)			
				{uni.list[[i]]<-tmp$unimesh
				#warplist[[i]]<-tmp$warp.mesh
				}
			if (is.character(rawout))
				{rawmesh<-rotmesh.onto(tmp$unimesh,proc$rotated[,,i],lm.data[,,i],scale=TRUE)$mesh
				mesh2ply(rawmesh,paste(rawout,lmnames[i],sep=""))
				}
			if (is.character(uniout))
				{mesh2ply(tmp$unimesh,paste(uniout,lmnames[i],sep=""))
				}
			}
		
		
		
		
	}
		if (keep)
		{names(uni.list)<-lmnames
		#names(proclist)<-lmnames
		}
		return(uni.list)
}

