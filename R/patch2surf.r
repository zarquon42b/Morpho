patch2surf<-function(patch,lm.ref,lm.tar,tarmesh,ray=FALSE,refmesh=NULL,tol=NULL,clean=TRUE)
	
{ 	if (is.null(tol))
		{tol<-1e3
		}
	est0<-tps3d(patch,lm.ref,lm.tar)
	if (ray)	
		{	if (is.null(refmesh))
				{stop("please enter reference mesh to estimate normals")
				}
			if (is.character(refmesh))
				{	refmesh_in<-file2mesh(refmesh)
				}
			est.mesh<-warp.mesh(refmesh_in,lm.ref,lm.tar,updateNormals=FALSE)
			mesh2ply(est.mesh,"est.mesh")
			patch_upd<-projRead(est0,"est.mesh.ply",clean=FALSE)
			system(paste("triray_project out_cloud.ply",tarmesh,tol,sep=" "))
			patchpro<-ply2mesh("out_ray.ply")
		
		}
	else
	{
	patchpro<-projRead(est0,tarmesh,readnormals=F)
	}
	if (clean)
		{unlink(c("est.mesh.ply","out_cloud.ply","out_ray.ply"))
		}
	return(patchpro)
}
