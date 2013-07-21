closemeshKD <- function(x,mesh,k=50,sign=FALSE,cores=1,method=0,...)
  {

    if(.Platform$OS.type == "windows")
       cores <- 1
    if (is.null(mesh$normals))
      {
        mesh <- adnormals(mesh)
      }
    if (is.matrix(x))
      {
        matr <- x
        x <- list()
        x$vb <- rbind(t(matr),1)
      }
    else
      {
        matr <- t(x$vb[1:3,])
      }
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(vb)[2]
    nit <- dim(it)[2]
    
    nmat<-dim(matr)[1]
     dif<-rep(0,nmat)
     fptr <- dif
    bary <- barycenter(mesh)
    clostInd <- mcNNindex(bary,matr,k=k,cores=cores,...)
    ##clostInd <- nn2(bary,matr,k=k,...)$nn.idx
    ##clostInd <- knnx.index(bary,matr,k=k,algorithm="kd_tree",...)
    storage.mode(k) <- "integer"
    clost <- c(0,0,0)
    storage.mode(fptr) <- "integer"
    storage.mode(it) <- "integer"
    storage.mode(nvb) <- "integer"    
    storage.mode(nit) <- "integer"
    storage.mode(method) <- "integer"
    storage.mode(nmat) <- "integer"
    storage.mode(matr) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    storage.mode(clostInd) <- "integer"
    outmatr <- matr
    region <- fptr
    out <- .Fortran("matr_meshKD",matr,nmat,vb,nvb,it,nit,clostInd,k,dif,fptr,outmatr,region,mesh$normals[1:3,],sign,t(matr),method)
    gc()
    x$vb[1:3,] <- t(out[[11]])
    x$quality <- out[[9]]
    x$normals <- rbind(out[[15]],1)
    x$ptr <- out[[10]]
    return(x)
  }
