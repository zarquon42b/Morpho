estShapeFull<-function(data,bone,skin,dowel,fixwarp,ref,n0,userot=NULL,shapelist=NULL,meshref=NULL,warp=FALSE,usebone=FALSE,real=FALSE,scale=FALSE)
{	
	if (is.null(userot))
		{userot<-bone	
		}
	n<-dim(data)[3]
	lm.est<-NULL
	meandev<-NULL
	all.proc<-data
	mesh<-NULL
	### perfom GPA on bone data and calculate weights###
	
	bone.proc<-mc.procGPA(data[bone,,],scale=scale)
	procrot<-bone.proc$rotated
	wcalc<-proc.weight(procrot,n0,ref)
	
	

	### perfom GPA on skin data and calculate preliminary estimate ###
	if (!usebone)	
	{
	softproc<-mc.procGPA(data[skin,,],scale=scale)$rotated
	
		### rotate all lm on superimposed skin config ###
		for (i in 1:n)
			{all.proc[,,i]<-rotonmat(data[,,i],data[skin,,i],softproc[,,i],scale=scale)
			}
		
	lms<-all.proc[,,wcalc$data$nr]
	lm.est<-matrix(0,dim(data)[1],3)
		if (i > 1)
			{	
			for (i in 1:n0)
				{lm.est<-lm.est+lms[,,i]*wcalc$data$weight[i]
				}
			}
		else 
			{lm.est<-lms
			}
		
		### rotate estimate back into original coordinate system according to userot information ###
	
		lm.rot<-rotonmat(lm.est,lm.est[userot,],data[userot,,ref],scale=scale)
		lm.real<-rotonmat(lm.est,lm.est[skin,],data[skin,,ref],scale=scale)
	}

	else
		{
		for (i in 1:n)
			{all.proc[,,i]<-rotonmat(data[,,i],data[bone,,i],procrot[,,i],scale=scale)
			}
		
		lms<-all.proc[,,wcalc$data$nr]
		lm.est<-matrix(0,dim(data)[1],3)
		if (i > 1)
		{	
		for (i in 1:n0)
			{lm.est<-lm.est+lms[,,i]*wcalc$data$weight[i]
			}
		}
	else 
		{lm.est<-lms
		}

		### rotate estimate back into original coordinate system according to userot information ###	

		lm.rot<-rotonmat(lm.est,lm.est[userot,],data[userot,,ref],scale=scale)
		lm.real<-rotonmat(lm.est,lm.est[skin,],data[skin,,ref],scale=TRUE)
	}
	
	
	
	
	### warp it on dowels while keeping landmarks defined in fixwarp ###
	lm.warp<-tps3d(lm.rot,lm.rot[c(fixwarp,dowel),],rbind(lm.rot[fixwarp,],data[dowel,,ref]))
	estim<-lm.warp[skin,]
	

	meandev<-mean(sqrt(diag(tcrossprod(lm.real[skin,]-data[skin,,ref]))))
	meandevreal<-mean(sqrt(diag(tcrossprod(lm.rot[skin,]-data[skin,,ref]))))
	diffs<-sqrt(diag(tcrossprod(lm.rot[skin,]-data[skin,,ref])))
	meandevreal<-mean(diffs)
	print(meandev)
	
	dimvb<-dim(shapelist[[1]]$vb)
	lm.est<-NULL
	vb<-matrix(0,3,dimvb[2])
	
	
		
	if (!is.null(shapelist))
		{
		dimvb<-dim(shapelist[[1]]$vb)
		
		refdat<-matrix(0,dim(meshref)[1],dim(meshref)[2])
		vb<-matrix(0,3,dimvb[2])
		for (i in 1:n0)
			{vb<-vb+shapelist[[wcalc$data$nr[i]]]$vb[1:3,]*wcalc$data$weight[i]
			refdat<-refdat+meshref[,,i]*wcalc$data$weight[i]
			}
			#print(names(shapelist)[wcalc$data$nr[i]])
			mesh<-shapelist[[1]]
			mesh$vb[1:3,]<-vb
			meshrot<-mesh
			#meshbone<-mesh
			if (real)
				{meshrot<-rotmesh.onto(mesh,refdat,data[skin,,ref],scale=TRUE)$mesh
				meshrot<-adnormals(meshrot)
				}
			else
				{meshrot<-rotmesh.onto(mesh,refdat,lm.rot[skin,],scale=TRUE)$mesh
				meshrot<-adnormals(meshrot)
				}
			if (warp)
			{			
			mesh<-mc.warp.mesh(mesh,refdat,estim)
			}
			else
				{mesh<-NULL
				}
			#meshbone<-mc.warp.mesh(meshrot,lm.rot[bone,],data[bone,,ref])
			#meshbon<-adnormals(mesh)	
			

			}	
		est.bone<-lm.rot[bone,]
		
	return(list(estim=estim,rot=lm.rot,weights=wcalc$data,meandev=meandev,meandevreal=meandevreal,mesh=mesh,meshrot=meshrot,real=lm.real[skin,],diffs=diffs,est.bone=est.bone))
	
	
}
