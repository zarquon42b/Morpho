unify.mesh<-function(mesh1,mesh2,data1,data2,ray=FALSE,tol=NULL,updateNormals=TRUE)
{ 	
	warp.m2<-warp.mesh(mesh1,data1,data2,updateNormals=TRUE)
	if (!ray)
		{uni.warp2<-mesh2mesh(warp.m2,mesh2)
		}
	else
		{updateNormals<-FALSE
		if (is.null(tol))
			{tol<-3*c.size(data1)/dim(data1)[1]}
                
		uni.warp2<-ray2mesh(warp.m2,mesh2,tol=tol)
		}
      	if(updateNormals)
		{cat("updating normals...\n")
		uni.warp2<-adnormals(uni.warp2)
		}
	return(list(unimesh=uni.warp2,warp.mesh=warp.m2))
}
