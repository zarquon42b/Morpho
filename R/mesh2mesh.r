mesh2mesh<-function(mesh1,tarmesh,clean=TRUE,cloud=FALSE,sign=FALSE)
{
  options <- NULL
  if (sign)
    {options <- paste(options,"--sign")
   }
  it <- NULL
  if (!is.character(tarmesh))
    {mesh2ply(tarmesh,"target")
     tarply <- "target.ply"
   }
  else
    {tarply <- tarmesh
   }
  if (!is.character(mesh1))
    {
      vert<-t(mesh1$vb[1:3,])
      projBack(vert,tarply)
      it <- mesh1$it
    }
  
  else
    {
      if (!cloud)
        {
          it <- file2mesh(mesh1)$it
        }
      system(paste("trimesh_project ",mesh1," ",tarply,options, sep=""))
    }
  outmesh <- ply2mesh("out_cloud.ply",readnormals=TRUE)
  class(outmesh) <- c("mesh3d","shapes3d")
  outmesh$vb <- rbind(outmesh$vb,1) 
  outmesh$it <- it
  if (!is.null(it))
    {
      outmesh <- adnormals(outmesh)
    }
  if (clean &&!is.character(tarmesh) )
    {
      if (!is.character(tarmesh))
        {unlink(c("target.ply"))
       }
      unlink("out_cloud.ply")
    }
  return(outmesh)
}
