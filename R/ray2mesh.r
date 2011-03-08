ray2mesh<-function(mesh1,tarmesh,tol=1,clean=TRUE,outname=NULL,readback=TRUE)
{	
	
	if (is.null(outname))
		{outname<-"project.mesh.ply"
		}
	mesh2ply(tarmesh,"target")
	mesh2ply(mesh1,"reference")
	system(paste("rayproject reference.ply target.ply ",tol,outname,sep=""))
	outmesh<-ply2mesh(outname,readnormals=TRUE)
	if (clean)
		{unlink(c("reference.ply","target.ply"))
		}
	if (readback)
		{return(outmesh)
		}
}
