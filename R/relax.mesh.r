relax.mesh <- function(mesh1,mesh2,ray=T,tol=1,split=1000,iter=1,lm=NULL,rhotol=0.7,sdmax=3,uselog=FALSE)
  {
    nlm <- NULL
    free <- NULL
    surface <- NULL
    lm1 <- lm1.norm <- NULL
    
    vb.lm <- NULL
    norm.lm <- NULL
    
    
    vb.m1 <- t(mesh1$vb[1:3,]) ### original vertices (of mesh1)
    norm.m1 <- t(mesh1$normals[1:3,]) ###original normals (of mesh1)
    
    if (ray) ## check if search along normals is selected
      {
        w.mesh<-ray2mesh(mesh1,mesh2,tol=tol)
        if (!is.null(lm))
          {lmini <- proj.read(lm,mesh1,readnormals=TRUE)  ### get normals from additional landmarks
          w.lm <- ray2mesh(lmini,mesh2,tol=tol)  ### throw landmarks on target along ray has to be worked on yet
         #  w.lm <- proj.read(t(lmini$vb[1:3,]),mesh2)
         }
      }
    else
      {
        w.mesh<-mesh2mesh(mesh1,mesh2)
        if (!is.null(lm))
          {lmini <- proj.read(lm,mesh1)  ### get normals from additional landmarks
           w.lm <- proj.read(lm,mesh2,readnormals=TRUE)  ###throw landmarks on target
         }
      }
### projected vertices and normal on target mesh ###

    vb1 <- t(w.mesh$vb[1:3,])
    norm1 <- t(w.mesh$normals[1:3,])
    dist1 <- w.mesh$quality
    k <- dim(vb1)[1]

    slideall <- NULL

    if (!is.null(lm))
      {
        nlm <- dim(lm)[1]
        lm1 <- t(lmini$vb[1:3,])### original landmarks on mesh1
        lm1.norm <- t(lmini$normals[1:3,]) ### original landmarks' normals (on mesh1)
        vb.lm <- t(w.lm$vb[1:3,]) ### projected landmarks
        norm.lm <- t(w.lm$normals[1:3,])### projected landmarks' normals
       
      }
   
   
   
        
### check for invalid normals of vertices ###
        
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
 ### check for invalid normals of landmarks ###
 if (!is.null(lm))
      {
        cnt <- 0
        degnlm <- NULL
        for ( i in 1:nlm)
      
          { if (prod(norm.lm[i,] == c(0,0,0))==1)
              {cnt <- cnt+1
               degnlm[cnt] <- i
             }
          }
        
        if (!is.null(degnlm))
          {
            lm1 <- lm1[-degnlm,]
            lm1.norm <- lm1.norm[-degnlm,]
            vb.lm<- vb.lm[-degnlm,]
            norm.lm <-  norm.lm[-degnlm,]

### update landmark number after removal of invalid normals
           
            nlm1 <- dim(lm1)[1]
          }      
      }
    

   
### make sure requested subsample size does not exceed number of vertices from the reference mesh ##

    if (k1 < split) 
      {
        split <- k1
      }
    
        sample.vb <- sample(1:k1) ### scramble vertex indices
        tmp <- sample.vb[(1:split)]#+(i*split)]
    
    dat <- rbind(vb.lm,vb1[tmp,])  ### subset of projected vertices
    #dist1 <- dist1[tmp] ### subset of vector of distances between projection and origin
    norm <- rbind(norm.lm,norm1[tmp,]) ### subset of normals of projected verticesp
    norm.m1 <- rbind(lm1.norm,norm.m1[tmp,]) ### subset of normals of original vertices
    vb.m1 <- rbind(lm1,vb.m1[tmp,]) ### subset of original vertices
    L <- CreateL(vb.m1)
    
### check if target vertices differ in their normal according to the value of rhotol
    rho <- NULL
    for (i in 1:split)
      {
        rho[i] <- angle.calc(norm[i,],norm.m1[i,])$rho
            
      }

### end angle check ####################################################### 


#### find projected points with high bending energy #############
 # print (dim(L$Lsubk))
  #  print(dim(dat))
    coeff <- L$Lsubk%*%dat
    coeff <- rbind(coeff,matrix(0,3,3))
    dif.be <- fx(vb.m1,dat,coeff)
    dif.be <- dif.be^2
    dif.be <- apply(dif.be,1,sum)
   
    
#### end bending energy check ############################     

    rhoex <- which(rho > rhotol)
    if (uselog)
      {
        logdif <- log(dif.be)
        diflogsd <-  mean(logdif)+sdmax*(sd(logdif))
        difex <- which(log(dif.be) > diflogsd)
      }
    else
      {
        difsd <- mean(dif.be)+sdmax*(sd(dif.be))
        difex <- which(dif.be > difsd)
      }
   
    tmp <- which(difex %in% rhoex)
    difex <- difex[-tmp]
    rhoex <- c(rhoex,difex)
   
    if (length(rhoex) > 0)
      {
        free <- rhoex
        surface <- (1:split)[-rhoex]
        }
    else
      {surface <- 1:split
     }
    gc()
#### start relaxation against reference ########
    
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

##### end relaxation #########################
    
    slideall <- dat
   
       
          
            mesh1.lm <- vb.m1

        w.mesh <- unify.mesh(mesh1,mesh2,mesh1.lm,slideall,ray=F,tol=tol)

    gc()
        return(list(norm=norm,dat=dat,mesh=w.mesh,slideall=slideall,dist1=dist1,rho=rho,dataslide=dataslide,rhoex=rhoex,dif.be=dif.be))
  }
    
