close <- function(point,matr,normals)
  {
    clost <- c(0,0,0)
    storage.mode(point) <- "double"
    storage.mode(matr) <- "double"
    storage.mode(normals) <- "double"
    storage.mode(clost) <- "double"
    
    a <- .Fortran("closest",point,matr,normals,clost)
    return(a)
  }

closemesh <- function(point,mesh)
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it
    fn <- dim(it)[2]
    clost <- c(0,0,0)
    dist_old<- 1e10
    distl <- 0
    for (i in 1:fn)
      {
       
        tmp <- close(point,t(vb[,it[,i]]),diag(1:3))[[4]]
        dist<-sqrt(sum((tmp-point)^2))
        distl[i] <- dist
#spheres3d(tmp,radius=0.3)
        
        if (dist < dist_old)
          {clost <- tmp
           dist_old<-dist
         }
      }
    return(list(clost=clost,dist=distl))
  }

closemeshf <- function(point,mesh)
  {
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(it)[2]
    nit <- dim(vb)[2]
    dif<-0
    
    clost <- c(0,0,0)
    normals <- diag(1:3)
    storage.mode(it) <- "integer"
        storage.mode(nvb) <- "integer"
    
    storage.mode(nit) <- "integer"

    storage.mode(normals) <- "double"
    storage.mode(point) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    out <- .Fortran("closemesh",point,vb,nvb,it,nit,normals,clost,dif)
    gc()
    return(out)
  }
closemeshf1 <- function(point,mesh)
  {
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(it)[2]
    nit <- dim(vb)[2]
    dif<-0
    fptr <- 0
    clost <- c(0,0,0)
    storage.mode(fptr) <- "integer"
    storage.mode(it) <- "integer"
    storage.mode(nvb) <- "integer"
    
    storage.mode(nit) <- "integer"

    
    storage.mode(point) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    out <- .Fortran("pt_mesh",point,vb,nvb,it,nit,clost,dif,fptr)
    gc()
    return(out)
  }
                
upsearch<-function(point,mesh)
  {
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(it)[2]
    nit <- dim(vb)[2]
    DAT<-matrix(0,nit,12)
      storage.mode(DAT)<-"double"

    storage.mode(it) <- "integer"
    storage.mode(nvb) <- "integer"
    storage.mode(nit) <- "integer"
    storage.mode(vb) <- "double"              
   a<- .Fortran("updateSearch",vb,nvb,it,nit,DAT)
    return(a)
  }
pt_tri<-function(point,mesh)
  {
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(it)[2]
    nit <- dim(vb)[2]
   }
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
mc.closemat <- function(matr,mesh,cores=getOption("cores")
  {
    if cores=NULL
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
    mcfun(i) <-function(i)
      {out <- .Fortran("matr_mesh",matr,nmat,vb,nvb,it,nit,dif,fptr,outmatr,region)[[9]]
       return(out)
     }
    mclist <- mclapply(mclist,closemat,mesh=mesh)
    gc()
    outmat <- NULL
    for (i in 1:cores)
      {outmat <- rbind(outmat,mclist[[i]])
                     }
    return(out)
  }
mesh_mesh<- function(mesh1,mesh,rhotol)9
    matr <- t(mesh1$vb[1:3,])
    if (is.null(mesh1$normals))
      {mesh1 <- adnormals(mesh1)
     }
    normals <- mesh1$normals
    tarnorm <- mesh$normals
    vb <- (mesh$vb[1:3,])
    it <- (mesh$it)
    nvb <- dim(it)[2]
    nit <- dim(vb)[2]
    
    nmat<-dim(matr)[1]
     dif<-rep(0,nmat)
     fptr <- 0
    clost <- c(0,0,0)
    storage.mode(fptr) <- "integer"
   
    storage.mode(it) <- "integer"
    storage.mode(nvb) <- "integer"
    storage.mode(normals) <- "double"
    storage.mode(tarnorm) <- "double"
    storage.mode(rhotol) <- "double"
    storage.mode(nit) <- "integer"
    storage.mode(nmat) <- "integer"
    storage.mode(matr) <- "double"
    storage.mode(point) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"

    out <- .Fortran("rho_mesh",matr,nmat,normals,vb,nvb,it,nit,tarnorm,dif,fptr)
    gc()
    return(out)
  }
pt_upmesh <- function(x,mesh)
{
  vb <- (mesh$vb[1:3,])
  it <- (mesh$it)
  nit <- dim(it)[2]
  nvb <- dim(vb)[2]
  dif<-0
  fptr <- 0
  clost <- c(0,0,0)
  storage.mode(fptr) <- "integer"
  storage.mode(it) <- "integer"
  storage.mode(nvb) <- "integer"
  storage.mode(nit) <- "integer"

    
    storage.mode(point) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    out <- .Fortran("pt_upmesh",point,vb,nvb,it,nit,clost,dif,fptr)
    gc()
    return(out)
}
