file2ply2 <- function(filein,fileout,curvature=TRUE,neighbour=2)
  {
    tmpmesh <- file2mesh(filein)
    mesh2ply2(tmpmesh,filename=fileout,curvature=curvature,neighbour=neighbour)
  }
