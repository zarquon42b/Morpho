barycenter <- function(mesh)
  {
    vb <- mesh$vb[1:3,]
    nvb <- dim(vb)[2]
    it <- mesh$it
    storage.mode(it) <- "integer"
    nit <- dim(it)[2]
    bary <- matrix(0,nit,3)
    storage.mode(bary) <- "double"

    out <- .Fortran("barycenter",vb,nvb,it,nit,bary)[[5]]
    return(out)
  }
