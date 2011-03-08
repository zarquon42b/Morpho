closemat <- function(matr,mesh)
  {
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(vb)[2]
    nit <- dim(it)[2]
   
    nmat<-dim(matr)[1]
     dif<-rep(0,nmat)
     fptr <- dif
    clost <- c(0,0,0)
    storage.mode(fptr) <- "integer"
   
    storage.mode(it) <- "integer"
    storage.mode(nvb) <- "integer"
    
    storage.mode(nit) <- "integer"
    storage.mode(nmat) <- "integer"
    storage.mode(matr) <- "double"
    storage.mode(point) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    outmatr <- matr
    region <- fptr
    out <- .Fortran("matr_mesh",matr,nmat,vb,nvb,it,nit,dif,fptr,outmatr,region)
    gc()
    return(out)
  }
