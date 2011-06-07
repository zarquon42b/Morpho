file2ply2 <- function(filein,fileout)
  {
    tmpmesh <- file2mesh(filein)
    mesh2ply2(tmpmesh,filename=fileout)
  }
