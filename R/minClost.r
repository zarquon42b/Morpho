minClost <- function(x,y)
  {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    dimX <- dim(x)
    dimY <- dim(y)
    out <- 1:dimX[1]

    dists <- .Fortran("minClost",x,dimX[1],y,dimY[1],out)
    return(dists[[5]])
  }
