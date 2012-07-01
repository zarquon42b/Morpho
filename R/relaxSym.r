relaxSym <- function(mesh1,mesh2,sigma,gamma=2,n,f)
{
  sigma0 <- sigma
  nmat=dim(mesh1$vb)[2]
  print(nmat)
  storage.mode(nmat) <- "integer"
  for (i in 1:n)
    {
      S0 <- t(mesh1$vb[1:3,])
      sigma <- (sigma0*f^(-i))^2
      S <- t(closemeshKD(mesh1,mesh2,sign=F)$vb[1:3,])
      M <-  t(closemeshKD(mesh2,mesh1,sign=F)$vb[1:3,])
      D1 <- t(mesh1$vb[1:3,]) - S
      D2 <- t(mesh2$vb[1:3,]) - M
      storage.mode(i) <- "integer"
      print(dim(S))
      time <- system.time(out <- .Fortran("displace_mesh_gauss",S,nmat,S,nrow(S),M,nrow(M),D1,D2,sigma,gamma,S*0));print(time)
      mesh1$vb[1:3,] <- t(S+out[[11]])
    }
  
  return(list(out=out,S=S,meshout=mesh1))
}
relaxSymtest <- function(mesh1,mesh2,sigma,gamma=2,n,f,oneway=F,k=1)
{
  sigma0 <- sigma
  M0 <- t(mesh2$vb[1:3,])
  S0 <- t(mesh1$vb[1:3,])
  sigma <- (sigma0*f^(-k))^2
  S <- t(closemeshKD(mesh1,mesh2,sign=F)$vb[1:3,])
  M <-  t(closemeshKD(mesh2,mesh1,sign=F)$vb[1:3,])
  D1 <- S-S0
  D2 <- M-M0
  storage.mode(S0) <- "double"
  storage.mode(M) <- "double"
  storage.mode(D1) <- "double"
  storage.mode(D2) <- "double"
  storage.mode(n) <- "integer"
                                        #out <- .Fortran("relax_pt",S[1,],S,nrow(S),M,nrow(M),D1,D2,sigma,gamma,k,c(0,0,0))
      time <- system.time(out <- .Fortran("displace_mesh_gauss",S0[n,],length(n),S0,nrow(S0),M,nrow(M),D1,D2,sigma,gamma,S0[n,]*0,oneway));print(time)
    
  
  return(list(out=out,S=S,addit=S0[n,]+out[[11]],M=M,S=S,D1=D1,D2=D2))
}
#require(Morpho)
#dyn.load("../Morpho/src/test.so")
