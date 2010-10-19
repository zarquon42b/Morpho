unify.mesh<-function(mesh1,mesh2,data1,data2,updateNormals=TRUE)
{ 	
	warp.m2<-warp.mesh(mesh1,data1,data2,updateNormals=FALSE)
	uni.warp2<-mesh2mesh(warp.m2,mesh2)
	
      	if(updateNormals)
		{cat("updating normals...\n")
		uni.warp2<-adnormals(uni.warp2)
		}
	return(list(unimesh=uni.warp2,warp.mesh=warp.m2))
}
