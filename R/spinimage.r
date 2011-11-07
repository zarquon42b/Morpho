spinimage <- function(mesh,O,bs,tol=NULL,rho=pi/2,iwidth=100)
  {
    ### tol=max tolerance for beta
    angcomp <- NULL
    VB <- mesh$vb[1:3,]
    normals <- mesh$normals[1:3,]
    dimvb <- dim(VB)[2]
    comp <- normals[1:3,O]
    if (!is.null(rho))
      {
        angtest <- function(i)
          {
            out <- angle.calc(normals[,i],comp)$rho
            return(out)
          }
        angcomp <- unlist(lapply(as.list(1:dimvb),angtest))
        rm <- which(abs(angcomp) < rho)
        VB <- VB[,rm]
      }
    
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
    if (is.null(tol))
      {tol <- mdist+1
     }
    else
      {
        bsmall <- which(abs(S0[,2]) <= tol)
        S0 <- S0[bsmall,]
      }
    bmax <- ceiling(2*tol)
    amax <- ceiling(mdist)
    imax <- iwidth
    #imax <- (2*bmax/bs)+1
    jmax <- imax
    #jmax <- (amax/bs)+1
    i <- floor((bmax-S0[,2])/bs)+1
    j <- floor(S0[,1]/bs)+1
   
    ij <- cbind(i,j)
     cleanij <- unique(c(which(i > iwidth),which(j > iwidth)))
    if (length(cleanij > 0))
      {ij <- ij[-cleanij,]
       S0 <- S0[-cleanij,]
     }
       ij <- ij[,]
   
   
    
    a <- S0[,1]-ij[,1]*bs
    b <- S0[,2]+ij[,2]*bs
                                        #bmax/2-ij[,2]*bs
   
    ab <- cbind(a,b)
    #ab <- ab[,]
    Sp <- matrix(0,imax,jmax)
    storage.mode(ij) <- "integer"
    storage.mode(imax) <- "integer"
    storage.mode(jmax) <- "integer"
    storage.mode(Sp) <- "double"
    storage.mode(bs) <- "double"
    storage.mode(ab) <- "double"    
    print(dim(ab))
    print(range(ij))
    print(dim(ij))
    print(dim(Sp))
    out <- NULL
    out <- .Fortran("spinimage",Sp,ij,dim(ij)[1],ab,imax,jmax,bs)[[1]]
    return(list(Sp=out,S0=S0,rm=rm,ij=ij))
  }
