icp <- function(mesh1,mesh2,iterations=3,scale=T,lm1=NULL,lm2=NULL,uprange=0.9,rhotol=pi,k=50)
  {

    if (!is.null(lm1))## perform initial rough registration
      {
        mesh1 <- rotmesh.onto(mesh1,lm1,lm2,scale=T)$mesh
        
      }
    
    
    for( i in 1:iterations)
      {
        proMesh <- closemeshKD(mesh1,mesh2,k=k) ## project mesh1 onto mesh2
        x1 <- vert2points(mesh1)
        x2 <- vert2points(proMesh)
        x2norm <- proMesh$normals
        x1norm <- mesh1$normals
        dists <- abs(proMesh$quality)
        
        ## check if normals angles are below rhotol
        normcheck <- rep(0,dim(x2norm)[2])
        storage.mode(normcheck) <- "double"
        storage.mode(x1norm) <- "double"
        storage.mode(x2norm) <- "double"
        ln <- dim(x2norm)[[2]]
        storage.mode(ln) <- "integer"
        normcheck <- .Fortran("angcheck",x1norm,ln,x2norm,normcheck)[[4]]
        goodnorm <- which(normcheck < rhotol)
        x1 <- x1[goodnorm,]
        x2 <- x2[goodnorm,]
        dists <- dists[goodnorm]

        ## check distances of remaining points and select points
        qud <- quantile(dists,probs=uprange)
        good <- which(dists < qud)
        mesh1 <- rotmesh.onto(mesh1,x1[good,],x2[good,],scale=scale)$mesh
        
      }
    return(mesh1)
  }
        

normcheck <- function(mesh1,mesh2)
  {
    
     x1norm <- mesh1$normals
     x2norm <- mesh2$normals
     ln <- dim(x2norm)[[2]]
     normcheck <- rep(0,dim(x2norm)[2])
     storage.mode(normcheck) <- "double"
     storage.mode(x1norm) <- "double"
     storage.mode(x2norm) <- "double"
     storage.mode(ln) <- "integer"
     normcheck <- .Fortran("angcheck",x1norm,ln,x2norm,normcheck)[[4]]
     return(normcheck)
   }
