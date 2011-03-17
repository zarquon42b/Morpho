place.patch <- function(dat.array,path,atlas.mesh,atlas.lm,lm,which.fix,prefix="skull_",tol=tol,split=1000,ray=T,outlines=NULL,SMvector=NULL,deselect=TRUE)
  {
    n <- dim(dat.array)[3]
    lm.dim <- dim(lm)
    out <- array(NA,dim=c(lm.dim,n))
    dimnames(out) <- dimnames(dat.array)
    
    L <- CreateL(atlas.lm)
        
    name <- dimnames(dat.array)[[3]]
    for(i in 1:n){

      
      tmp.name <- paste(path,prefix,name[i],".ply",sep="")
      tmp.data <- proj.read(dat.array[,,i],tmp.name,readnormals=TRUE)
      
      U<-calcTang_U_s(t(tmp.data$vb[1:3,]),t(tmp.data$normals[1:3,]),SMvector=SMvector,outlines=outlines,surface=NULL,deselect=deselect)
      slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
      tps.lm <- tps3d(lm,atlas.lm,slide)
      tps.lm <- proj.read(tps.lm,tmp.name,readnormals=FALSE)
      
      dat <- tps.lm
     
      out[,,i] <- dat

    }
    return(out)
  }
      
