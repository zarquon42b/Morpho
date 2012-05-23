bone2skin<-function(data,bonemesh,skinmesh)
{
  #out<-list()
  #out$vb<-t(data)
  clB <- clS <- FALSE
  if (!is.character(bonemesh))
    {
      mesh2ply(bonemesh,file="bonemesh")
      bonemesh <- "bonemesh.ply"
      clB=TRUE
    }
  if (!is.character(skinmesh))
    {
      mesh2ply(skinmesh,file="skinmesh")
      skinmesh <- "skinmesh.ply"
      clS=TRUE
    }
  projBack(data,bonemesh)
  system(paste("rayproject out_cloud.ply ", skinmesh," -t 100",sep=""))
  out<-ply2mesh("project.mesh.ply",readnormals=TRUE)
  unlink(c("project.mesh.ply","out_cloud.ply"))

  ##clean up
  if (clB)
    {unlink( "bonemesh.ply")
    }
  if (clS)
    {unlink( "skinmesh.ply")
   }
  return(out)
}
