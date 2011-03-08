relax.mesh <- function(mesh1,mesh2,ray=T,tol=1,split=1000,iter=1,lm=NULL,rhotol=0.7)
  {
    nlm <- NULL
    free <- NULL
    surface <- NULL
    vb.m1 <- t(mesh1$vb[1:3,])
    norm.m1 <- t(mesh1$normals[1:3,])
    
    if (ray) ## check if search along normals is selected
      {
        w.mesh<-ray2mesh(mesh1,mesh2,tol=tol)
        if (!is.null(lm))
          {lmini <- proj.read(lm,mesh1,readnormals=TRUE)
           w.lm <- ray2mesh(lmini,mesh2,tol=tol)
         }
      }
    else
      {
        w.mesh<-mesh2mesh(mesh1,mesh2)
        if (!is.null(lm))
          {lmini <- proj.read(lm,mesh1)
           w.lm <- proj.read(lmini,mesh2,readnormals=TRUE)
         }
      }
### projected vertices and normal on target mesh ###

    vb1 <- t(w.mesh$vb[1:3,])
    norm1 <- t(w.mesh$normals[1:3,])
    dist1 <- w.mesh$quality
      
    if (!is.null(lm))
      {
        nlm <- dim(lm)[1]
        vb.m1 <- rbind(lm,vb.m1)
        
        vb1 <- rbind(t(w.lm$vb[1:3,],vb1))
        norm1 <- rbind(t(w.lm$normals[1:3,],norm1))
        
      }
   
    k <- dim(vb1)[1]
    slideall <- NULL
        
### check for invalid normals ###
        
    cnt <- 0
    degnorm <- NULL
    for ( i in 1:k)
      
      { if (prod(norm1[i,] == c(0,0,0))==1)
         {cnt <- cnt+1
          degnorm[cnt] <- i
        }
     }
  
    if (!is.null(degnorm))
      {
        vb.m1 <- vb.m1[-degnorm,]
        norm.m1 <-  norm.m1[-degnorm,]
        vb1 <- vb1[-degnorm,]
        norm1 <- norm1[-degnorm,]
        dist1 <- dist1[-degnorm]
      }
### update vertex number after removal of invalid normals
        k1 <- dim(vb1)[1]
### make sure requested subsample size does not exceed number of vertices from the reference mesh ##
    if (k1 < split) 
      {
        split <- k1
      }
        sample.vb <- sample(1:k1)
        
        tmp <- sample.vb[(1:split)]#+(i*split)]
    if (!is.null(lm))
      {
        tmp <- c(1:nlm,(tmp+nlm))
      }
    dat <- vb1[tmp,]
    dist1 <- dist1[tmp]
    norm <- norm1[tmp,]
    norm.m1 <- norm.m1[tmp,]
    L <- CreateL(vb.m1[tmp,])
### check if target vertices differ in their normal according to the value of rhotol
    rho <- NULL
    for (i in 1:split)
      {
        rho[i] <- angle.calc(norm[i,],norm.m1[i,])$rho
            
      }
        #print((rho))
    rhoex <- which(rho > rhotol)
   
    if (length(rhoex) > 0)
      {
        free <- rhoex
        surface <- (1:split)[-rhoex]
        }
    else
      {surface <- 1:split
     }
    gc()
        for (i in 1:iter)
          { 
            cat(paste("iteration",i,"\n"))
            U<-calcTang_U_s(dat,norm,SMvector=1:split,surface=surface,free=free)
            gc()
            dataslide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
            gc()
            pro <- proj.read(dataslide,mesh2)
            dat <- t(pro$vb)
            norm <- t(pro$normals)
  
            gc()
            if (iter > 1)
              {
                rho <- NULL
                for (i in 1:split)
                  {
                    rho[i] <- angle.calc(norm[i,],norm.m1[i,])$rho
                    
                  }
                                        #print((rho))
                rhoex <- which(rho > rhotol)
                if (length(rhoex > 0))
                  {
                    free <- rhoex
                    surface <- (1:split)[-rhoex]
                  }
              
                else
                  {surface <- 1:split
                 }
                gc()
              }
          }
        
        slideall <- dat
   
        if (!is.null(degnorm))
          { mesh1.lm <- t(mesh1$vb[1:3,-degnorm])[sample.vb[1:(1*split)],]
          }
        else
          {
            mesh1.lm <- t(mesh1$vb[1:3,])[sample.vb[1:(1*split)],]
          }
        w.mesh <- unify.mesh(mesh1,mesh2,mesh1.lm,slideall,ray=F,tol=tol)
  # w.mesh$vb[1:3,] <- t(dataslide)
 #  w.mesh <- adnormals(w.mesh)
    
        return(list(norm=norm,dat=dat,mesh=w.mesh,slideall=slideall,dist1=dist1,rho=rho,dataslide=dataslide,rhoex=rhoex))
  }
    
