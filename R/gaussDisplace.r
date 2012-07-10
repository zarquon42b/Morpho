gaussDisplace <- function(mesh1,mesh2,sigma,gamma=2,W0,f,oneway=F,k=1,nh=NULL,tol=0,cores=detectCores(),k0=50,...)
{
  t0 <- Sys.time()
  sigma0 <- sigma
  M0 <- t(mesh2$vb[1:3,])
  S0 <- t(mesh1$vb[1:3,])
  sigma <- (sigma0*f^(-k))^2
  S <- vert2points(closemeshKD(mesh1,mesh2,sign=F,k=k0))
  t1 <- Sys.time()
    
   if (oneway)
    {
      M <- vert2points(mesh2)
    }
  else
    {
      #M <-  projRead(vert2points(mesh2),mesh1,sign=F,smooth=F,ignore.stdout=F,readnormals=F)
      M <- vert2points(closemeshKD(mesh2,mesh1,sign=F,k=k0))
    }
    t2 <- Sys.time()
  if (!is.null (nh))
    {
      #if (nh < 20)
       # {
        #  clostIndW <- nn2(S,W0,k=nh,searchtype="priority",...)$nn.idx
         # clostIndP <- nn2(M,W0,k=nh,searchtype="priority",...)$nn.idx
       # }
      #else
       # {
          clostIndW <- mcNNindex(S,W0,k=nh,...)
          clostIndP <- mcNNindex(M,W0,k=nh,...)
        #}
    }
  t3 <- Sys.time()
  D1 <- S-S0
  D2 <- M-M0
  storage.mode(clostIndW) <- "integer"
  storage.mode(clostIndP) <- "integer"
  #print(t1-t0)

  #print(t2-t1)
  #print(t3-t2)
 
  storage.mode(S0) <- "double"
  storage.mode(M) <- "double"
  storage.mode(D1) <- "double"
  storage.mode(D2) <- "double"
  storage.mode(nh) <- "integer"
  tol <- tol^2
  time <- system.time(out <- .Fortran("displace_mesh_gauss",W0,nrow(W0),S0,nrow(S0),M,nrow(M),D1,D2,sigma,gamma,oneway,clostIndW,nh,clostIndP,tol=tol))
  #print(time)
  return(list(out=out,S=S,addit=W0+out[[1]],M=M,S=S,D1=D1,D2=D2))
}

gaussDisplMesh3d <- function(mesh1,mesh2,iterations=10,smooth=NULL,smoothit=10,sigma=20,gamma=2,f=1.2,oneway=F,k=1,tol=0,lm1=NULL,lm2=NULL,icp=FALSE,icpiter=3,uprange=0.95,rhotol=1,nh=50,toldist=0,cores=detectCores(),k0=50,patch=NULL,repro=FALSE,...)
  {
    if (!is.null(patch))
      { ## append landmarks to meshes vertices
        mesh1$vb <- cbind(mesh1$vb,rbind(t(patch),1))
        colnames(mesh1$vb) <- NULL
        mdim <- dim(mesh1$vb)
        cols <- c((mdim[2]+1-dim(patch)[1]):mdim[2])
      }
    if (icp)
      {
        cat("performing icp matching\n")
        mesh1 <- icp(mesh1,mesh2,lm1=lm1,lm2=lm2,uprange=uprange,rhotol=rhotol,iterations=icpiter)
      }
    for (i in 1:iterations)
      {
        time0 <- Sys.time()
        
        if (!is.null(smooth))
          { if (i %% smooth == 0)
              {
                cat("smoothing step\n")
                mesh2ply(mesh1,"dump01")
                command <- paste("trismooth dump01.ply -it ",smoothit," dump01.ply",sep="") 
                system(command)
                mesh1 <- ply2mesh("dump01.ply")
                unlink("dump01.ply")
              }
            tmp <- gaussDisplace(mesh1,mesh2,sigma=sigma,gamma=gamma,f=f,W0=vert2points(mesh1),nh=nh,k=i,tol=toldist,cores=cores,k0=k0)
            mesh1$vb[1:3,] <- t(tmp$addit)
            if (!is.null(patch) && repro)
              {
                mesh1$vb[1:3,cols] <- t(closematKD(t(mesh1$vb[1:3,cols]),mesh1)[[11]])
              }
          }
        time1 <- Sys.time()
        cat(paste("completed iteration",i, "in", time1-time0, "seconds\n"))
        cat("\n")
      }
    if (!is.null(patch))
      invisible(list(mesh=mesh1,patch=vert2points(mesh1)[cols,]))
    else
      invisible(mesh1)
  }

    
