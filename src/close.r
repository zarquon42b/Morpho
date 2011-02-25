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
    vb <- t(mesh$vb[1:3,])
    it <- t(mesh$it)
    nvb <- dim(it)[1]
    nit <- dim(vb)[1]
    dif<-0
    clost <- c(0,0,0)
    normals <- diag(1:3)
    storage.mode(it) <- "integer"
    storage.mode(normals) <- "double"
    storage.mode(point) <- "double"
    storage.mode(vb) <- "double"
    storage.mode(dif) <- "double"
    storage.mode(clost) <- "double"
    out <- .Fortran("closemesh",point,vb,nvb,it,nit,normals,clost,dif)
   
    return(out)
  }
                

                   
                    
