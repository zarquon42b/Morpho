ang <- function(x,y)
  {
    storage.mode(x)<-"double"
    storage.mode(y)<-"double"
    lx <- length(x)
    ly <- length(y)
          storage.mode(lx)<-"integer"
    storage.mode(ly)<-"integer"
    rho <- 0
    storage.mode(rho)<-"double"

    a <- .Fortran("angcal",x,lx,y,ly,rho)
    return(a)
  }
