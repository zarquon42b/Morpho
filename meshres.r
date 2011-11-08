meshres <- function(mesh) UseMethod("meshres")
meshres.mesh3d <- function(mesh)
  {
    res <- 0
    VB <- mesh$vb[1:3,]
    nvb <- dim(VB)[2]
    IT <- mesh$it[1:3,]
    nit <- dim(IT)[2]
    storage.mode(VB) <- "double"
    storage.mode(IT) <- "integer"

    res <- .Fortran("meshres",VB,nvb,IT,nit,res)[[5]]
    return(res)
  }
    
