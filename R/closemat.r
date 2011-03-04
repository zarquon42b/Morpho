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
mc.closemat <- function(matr,mesh,cores=getOption("cores"))
  {
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
    storage.mode(point) <- "double"
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


mesh_mat <- function(mesh1,mesh2,toldist=1,rhotol=0.5)
{
   if (is.null(mesh1$normals))
     {mesh1 <- adnormals(mesh1)
     }
   if (is.null(mesh1$normals))
     {mesh2 <- adnormals(mesh2)
     }
   matr <- t(mesh1$vb[1:3,])
   matnorm <- t(mesh1$normals)
   VB <- (mesh2$vb[1:3,])
   IT <- (mesh2$it)
   nvb <- dim(VB)[2]
   tarnorm <- mesh2$normals
   nit <- dim(IT)[2]  
   nmat<-dim(matr)[1]
   dif<-rep(0,nmat)
   fptr <- dif
   clost <- c(0,0,0)
   storage.mode(rhotol) <- "double"
   storage.mode(toldist) <- "double"
   storage.mode(matr) <- "double"
   storage.mode(matnorm) <- "double"
   storage.mode(tarnorm) <- "double"
   storage.mode(dif) <- "double"
   storage.mode(clost) <- "double"
   storage.mode(VB) <- "double"
   storage.mode(fptr) <- "integer"
    storage.mode(IT) <- "integer"
    storage.mode(nvb) <- "integer" 
    storage.mode(nit) <- "integer"
    storage.mode(nmat) <- "integer"
   outmatr <- matr
   region <- dif
   #out <- .Fortran("mesh_mesh",matr,matnorm,nmat,VB,tarnorm,nvb,IT,nit,dif,fptr,outmatr)
   out <- .Fortran("mesh_mesh",matr,nmat,VB,nvb,IT,nit,dif,fptr,outmatr,region,matnorm,tarnorm,rhotol,toldist)
    gc()
   outmesh <- mesh1
   outmesh$vb[1:3,] <-t(out[[9]])
   outmesh <- adnormals(outmesh)
   return(list(out=out,mesh=outmesh))
 }
   
   
