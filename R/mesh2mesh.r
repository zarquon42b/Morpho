mesh2mesh<-function(mesh1,tarmesh,clean=TRUE)
{	mesh2ply(tarmesh,"target")
	vert<-t(mesh1$vb[1:3,])
	proj.back(vert,"target.ply")
	imp.vert<-ply2mesh("out_cloud.ply")
	outmesh<-mesh1
	outmesh$vb[1:3,]<-t(imp.vert)
	
	if (clean)
		{unlink(c("out_cloud.ply","target.ply"))
		}
	return(outmesh)
}
