warp.mesh<-function(mesh,matr,matt,lambda=0,updateNormals=TRUE)
{
      t00<-Sys.time()	
      vert<-t(mesh$vb[1:3,])
      
      cat("calculating spline...\n")
      warp<-tps3d(vert,matr,matt,lambda=lambda)
      mesh$vb<-rbind(t(warp),1)
      mesh$normals<-NULL
      
      
      testref<-rotonto(matr,matt)$reflect
      
  	if(testref == 1)
  		
      {mesh<-conv2backf(mesh)
        #cat("reflection is involved\n")
	}
      
      if(updateNormals)
	{cat("updating normals...\n")
	mesh<-adnormals(mesh)
	
	}
      t10<-Sys.time()
	#cat(paste("warping finished in", format(difftime(t10,t00,units="mins")[[1]]," minutes",sep="")))
      return(mesh)
}
      
