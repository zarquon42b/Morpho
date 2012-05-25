spinimage <- function(mesh,O,bs,rho=pi/2,resA=NULL,resB=NULL)
  {
   
    ### tol=max tolerance for beta
    angcomp <- NULL
    VB <- mesh$vb[1:3,]
    normals <- mesh$normals[1:3,]
    dimvb <- dim(VB)[2]
    comp <- normals[1:3,O]
    if (!is.null(mesh$it))
      {
        reso <- meshres(mesh)
        bs <- bs*reso
      }
### check angles
    if (!is.null(rho))
      {
        angcomp <- rep(0,dimvb)
        mat <- t(mesh$normals[1:3,])
        storage.mode(angcomp) <- "double"
        storage.mode(mat) <- "double"
        dimat <- dim(mat)[1]
        dimat2 <- dim(mat)[2]
        angcomp <- .Fortran("angcomp",angcomp,comp,mat,dimat,dimat2)[[1]]
        rm <- which(abs(angcomp) < rho)
        VB <- VB[,rm]
      }
### end check angles  
    
    bbox <- expand.grid(range(VB))
    mdist <- max(dist(bbox))    
    dimvb <- dim(VB)[2]
    p <- mesh$vb[1:3,O]
    n <- mesh$normals[1:3,O]
    S0 <- matrix(0,dimvb,2)
    storage.mode(S0) <- "double"
    storage.mode(VB) <- "double"
    storage.mode(n) <- "double"
### create spinmap
    S0 <- .Fortran("spinmap",S0,VB,p,n,dimvb)[[1]]
    
### create spinimage
    if (is.null(resB))
      {
        W <- mdist
        resB <- ceiling((2*W/bs)+1)        
     }
    else
      {
        W <-((resB-2)*bs)
        bsmall <- which(abs(S0[,2]) <= W/2)
        S0 <- S0[bsmall,]
      }
    if (is.null(resA))
      {
        resA <- ceiling((mdist/bs)+1)
      }
    S0[,2] <- W/2-S0[,2]
    amax <- mdist
    imax <- resB  
    jmax <- resA
    i <- floor((S0[,2])/bs)
    j <- floor(S0[,1]/bs)
    jclean <- which(j < (resA-1))
    ij <- cbind(i,j)   
    ij <- ij[jclean,]
    a <- S0[jclean,2]-(ij[,1])*bs
    b <- S0[jclean,1]-(ij[,2])*bs                                      #bmax/2-ij[,2]*bs
    ab <- cbind(a,b)
   # print(range(ab[,1]))
   #     print(range(ab[,2]))

    ij <- ij+1
    #ab <- ab[,]
    Sp <- matrix(0,imax,jmax)
    storage.mode(ij) <- "integer"
    storage.mode(imax) <- "integer"   
    storage.mode(jmax) <- "integer"
    storage.mode(Sp) <- "double"
    storage.mode(bs) <- "double"
    storage.mode(ab) <- "double"    
   # print(dim(ab))
   # print(range(ij))
   # print(dim(ij))
   # print(dim(Sp))
    out <- NULL
    out <- .Fortran("spinimage",Sp,ij,dim(ij)[1],ab,imax,jmax,bs)[[1]]
    return(list(Sp=out,S0=S0,rm=rm,ij=ij,ab=ab))
  }
spinimageLM <- function(mat,mesh,bs,rho=pi/2,resA=NULL,resB=NULL)
  {
    i <- NULL
    registerDoParallel() ##register parallel backend
    pro <- projRead(mat,mesh,readnormals=TRUE)
    vn <- dim(mesh$vb)[2]
    nlm <- dim(mat)[1]
    out <- list()
    tmpmesh <- mesh
    tmpmesh$vb <- cbind(tmpmesh$vb,rbind(pro$vb,1))
    tmpmesh$normals <- cbind(tmpmesh$normals,pro$normals)
    mcspin <- function(x)
      {
        ret <- spinimage(tmpmesh,x,bs=bs,rho=rho,resA=resA,resB=resB)$Sp
        return(ret)
      }
        
        
    out <- foreach (i= ((1:nlm)+vn)) %dopar% mcspin(i)
      return(out)
  }
 spinimageMesh <- function(mesh,bs,rho=pi/2,resA=NULL,resB=NULL)
  {
     i <- NULL
    registerDoParallel() ##register parallel backend
    
    vn <- dim(mesh$vb)[2]
    mcspin <- function(x)
      {
        ret <- spinimage(mesh,x,bs=bs,rho=rho,resA=resA,resB=resB)$Sp
        return(ret)
      }
    out <- foreach (i= (1:vn)) %dopar% mcspin(i)
    return(out)
  }
