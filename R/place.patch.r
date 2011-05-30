place.patch <- function(dat.array,path,atlas.mesh,atlas.lm,patch,which.fix,prefix=NULL,tol=5,ray=T,outlines=NULL,SMvector=NULL,deselect=TRUE,inflate=NULL)
  {
    n <- dim(dat.array)[3]
    k <- dim(dat.array)[1]
    patch.dim <- dim(patch)[1]
    out <- array(NA,dim=c((patch.dim+k),3,n))
    dimnames(out)[[3]] <- dimnames(dat.array)[[3]]
    
    L <- CreateL(atlas.lm)
    L1 <- CreateL(rbind(atlas.lm,patch))
    name <- dimnames(dat.array)[[3]]

    for(i in 1:n)
      {
        tmp.name <- paste(path,prefix,name[i],".ply",sep="")
        tmp.data <- proj.read(dat.array[,,i],tmp.name,readnormals=TRUE)

### relax existing curves against atlas ###
        if (!is.null(outlines))
          {
            sm <- SMvector
            U<-calcTang_U_s(t(tmp.data$vb[1:3,]),t(tmp.data$normals[1:3,]),SMvector=SMvector,outlines=outlines,surface=NULL,deselect=deselect)
            slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
            tps.lm <- tps3d(patch,atlas.lm,slide)
          }
        else
          
          {
            sm <- 1:k
            slide <- t(tmp.data$vb[1:3,])           
            tps.lm <- tps3d(patch,atlas.lm,t(tmp.data$vb[1:3,]))
            

          }
### use for mullitlayer meshes to avoid projection inside
        if (!is.null(inflate))
          {
          atlas.warp <- warp.mesh(atlas.mesh,atlas.lm,slide)
          tps.lm <- proj.read(tps.lm,atlas.warp,readnormals=TRUE,smooth=FALSE)
          tps.lm$vb[1:3,] <- tps.lm$vb[1:3,]+inflate*tps.lm$normals[1:3,] ###inflate outward along normals
          tps.lm <- ray2mesh(tps.lm,tmp.name,inbound=TRUE,tol=tol) ### deflate in opposite direction
          relax <- rbind(slide,t(tps.lm$vb[1:3,]))
          normals <- rbind(slide,t(tps.lm$normals[1:3,]))
              
          U1 <-calcTang_U_s(relax,normals,SMvector=sm,outlines=outlines,surface=c((k+1):(patch.dim+k)),deselect=deselect)                       
          tps.lm <- calcGamma(U1$Gamma0,L1$Lsubk3,U1$U,dims=3)$Gamatrix[c((k+1):(patch.dim+k)),]
          tps.lm <- proj.read(tps.lm,tmp.name,readnormals=FALSE)
        }
### just project warped patch on surface (suitable for singlelayer meshes)
        else
          {
            tps.lm <- proj.read(tps.lm,tmp.name,readnormals=FALSE)
          }
        out[,,i] <-rbind(dat.array[,,i],tps.lm)
      }
    return(out)
  }

