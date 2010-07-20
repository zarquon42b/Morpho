warp.mesh<-function(mesh,matr,matt,updateNormals=TRUE)
{
      
      vert<-t(mesh$vb[1:3,])
      
      cat("calculating spline...\n")
      warp<-tps3d(vert,matr,matt)
      mesh$vb<-rbind(t(warp),1)
      mesh$normals<-NULL
      
      
      sv1<-svd(t(matt)%*%matr)
      if (sign(det(sv1$v))<0)
      {mesh<-conv2backf(mesh)
        cat("reflection is involved\n")}
      
      if(updateNormals)
	{cat("updating normals...\n")
	mesh<-adnormals(mesh)
	
	}
      
      return(mesh)
}
      
