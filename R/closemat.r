closemat <- function(matr,mesh,sign=FALSE)
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
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    outmatr <- matr
    region <- fptr
    out <- .Fortran("matr_mesh",matr,nmat,vb,nvb,it,nit,dif,fptr,outmatr,region,mesh$normals[1:3,],sign,t(matr))
    gc()
    out <- out[c(7:10,13)]
    names(out) <- c("dist","fptr","clost","region","clostnorm")
    return(out)
  }

closemeshKD <- function(x,mesh,k=50,sign=FALSE,cores=detectCores(),method=0,...)
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
mcClosemat <- function(matr,mesh,cores=detectCores())
  {
    if(.Platform$OS.type == "windows")
      cores <- 1
    if (is.null(cores))
      {cores=3
     }
    mclist <- list()
    
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(vb)[2]
    nit <- dim(it)[2]
    nmat<-dim(matr)[1]
    dif<-rep(0,nmat)
    fptr <- dif
    clost <- c(0,0,0)

    iter <- as.integer(nmat/cores)
    for (i in 1:(cores-1))
      {mclist[[i]] <- matr[(1:iter)+((i-1)*iter),]
     }
    mclist[[cores]] <-  matr[-c(1:((cores-1)*iter)),]
   

    # set variable storage mode
    storage.mode(fptr) <- "integer"
    storage.mode(it) <- "integer"
    storage.mode(nvb) <- "integer" 
    storage.mode(nit) <- "integer"
    storage.mode(nmat) <- "integer"
    storage.mode(matr) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"

    outmatr <- matr
    region <- fptr
    mcfun1 <-function(x)
      {out <- closemat(x,mesh)
       return(out)
     }
    mclist <- mclapply(mclist,mcfun1)
    gc()
    outmat <- NULL
    dist <- NULL
    for (i in 1:cores)
      {outmat <- rbind(outmat,mclist[[i]][[9]])
       dist <- c(dist,mclist[[i]][[7]])
     }
    return(list(clost=outmat,dist=dist))
  }
