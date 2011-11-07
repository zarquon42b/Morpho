spinmap <- function(mesh,O)
  {
    VB <- mesh$vb[1:3,]
    dimvb <- dim(mesh$vb[1:3,])[2]
    p <- mesh$vb[1:3,O]
    n <- mesh$normals[1:3,O]
    S0 <- matrix(0,dimvb,2)
    storage.mode(S0) <- "double"
    storage.mode(VB) <- "double"
    storage.mode(n) <- "double"
    out <- .Fortran("spinmap",S0,VB,p,n,dimvb)[[1]]
    return(out)
  }   
    
