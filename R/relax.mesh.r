relax.mesh <- function(mesh1,mesh2,ray=T,tol=1,split=1000,iter=1,lm=NULL)
  {
    
    
        vb.m1 <- t(mesh1$vb[1:3,])

        if (ray)
          {
            w.mesh<-ray2mesh(mesh1,mesh2,tol=tol)
            if (!is.null(lm))
              {lmini <- proj.read(lm,mesh1)
               w.lm <- ray2mesh(lmini,mesh2,tol=tol)
             }
          }
        else
          {
            w.mesh<-mesh2mesh(mesh1,mesh2)
            if (!is.null(lm))
              {lmini <- proj.read(lm,mesh1)
               w.lm <- proj.read(lmini,mesh2)
             }# w.mesh <- adnormals(w.mesh)
          }
        vb1 <- t(w.mesh$vb[1:3,])
        norm1 <- t(w.mesh$normals[1:3,])
      
    if (!is.null(lm))
      {
        vb.m1 <- rbind(lm,vb.m1)
       # lm2mesh <- proj.read(lm,mesh2)
        vb1 <- rbind(vb1,t(w.lm$vb[1:3,]))
        norm1 <- rbind(norm1,t(w.lm$normals[1:3,]))
        
      }
   
    k <- dim(vb1)[1]
    slideall <- NULL
    # check for invalid normals
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
        vb1 <- vb1[-degnorm,]
        norm1 <- norm1[-degnorm,]
      }
    ### update vertex number after removal of ivalid normals
    k1 <- dim(vb1)[1]
   
    sample.vb <- sample(1:k1)
    tmp <- sample.vb[(1:split)]#+(i*split)]
    dat <- vb1[tmp,]
    norm <- norm1[tmp,]
     L <- CreateL(vb.m1[tmp,])
    gc()
    for (i in 1:iter)
      { 
        cat(paste("iteration",i,"\n"))
        U<-calcTang_U_s(dat,norm,SMvector=1:split,surface=1:split)
        gc()
        dataslide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
        gc()
        pro <- proj.read(dataslide,mesh2)
        dat <- t(pro$vb)
        norm <- t(pro$normals)
  
        gc()
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
    
    return(list(norm=norm,dat=dat,mesh=w.mesh,slideall=slideall))
  }
    
