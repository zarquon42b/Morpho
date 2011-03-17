place.lm <- function(data,path,atlas.mesh,atlas.lm,lm,which.fix,prefix="skull_",tol=tol,split=1000,ray=T)
  {
    n <- dim(data)[3]
    lm.dim <- dim(lm)
    out <- array(NA,dim=c(lm.dim,n))
    dimnames(out) <- dimnames(data)
    
    name <- dimnames(data)[[3]]
    for(i in 1:n){
      tmp.name <- paste(path,prefix,name[i],".ply",sep="")
      tmp.mesh <- file2mesh(tmp.name)
      tmp <- relax.lm(atlas.mesh,atlas.lm,tmp.mesh,data[,,i],lm=lm,tol=tol,split=split,ray=ray)
      dat <- tmp$dat[1:lm.dim[1],]
      dat[which.fix,] <- data[,,i]
      out[,,i] <- dat

    }
    return(out)
  }
      
