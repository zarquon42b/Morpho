
map.mesh <- function(mesh1,lm1,mesh2,lm2,tol=1e-3,it=2,overlap=0.8,raytol=NULL,strict=TRUE,n=NULL,uselm=TRUE,subset=NULL)
  {
    round <- 0
    p <- 1e10
    rsme <- 1e11
   
    rot <- rotmesh.onto(mesh1,lm1,lm2)
    tmp.mesh1 <- rot$mesh
     open3d()
    wire3d(mesh1,col=3)
    if (interactive())
      { f <-  select3d()
        subset <- t(mesh1$vb[1:3,])
        subset <- subset[which(f(subset)),]
        print(subset)
      }
     if (!is.null(subset))
       {rotsub <- rotonmat(subset,lm1,lm2)
        tmp.mesh <- proj.read(subset,tmp.mesh1)#rotsub$mesh
        tmp.mesh$vb <- rbind(tmp.mesh$vb,1)
      }
    tmp.lm <- rot$yrot
    
    ref.lm <- t(tmp.mesh$vb[1:3,])
    lmdim <- dim(lm1)[1]
   
    while(p > tol && round < it )
      {
        round <- round+1
        rsme_old <- rsme
        tmp.mesh_old <- tmp.mesh
        tmptar <- rbind(tmp.lm,t(tmp.mesh$vb[1:3,]))

        if (is.null(raytol))
          {
            tmptar <- proj.read(tmptar,mesh2)
          }

        else

          {tmptar <- proj.read(tmptar,tmp.mesh)
           
            tmptar <-  ray2mesh(tmptar,mesh2,tol=raytol,strict=strict)
         }
        
        tar.lm <- t(tmptar$vb[1:3,])
        ref.lm <- rbind(tmp.lm,t(tmp.mesh$vb[1:3,]))
        tar.norm <- t(tmptar$normals[1:3,])
        quant <- quantile(tmptar$quality,probs=overlap)
        good <- which(tmptar$quality[-c(1:lmdim)] < quant)
        good <- good + lmdim
        tar.lm <- tar.lm[c(1:lmdim,good),]
        ref.lm <- ref.lm[c(1:lmdim,good),]
        tar.norm <-  tar.norm[c(1:lmdim,good),]

        print(1)
### use a subset of vertices and relax against reference
        if (!is.null(n)) 
          {
            if (n > dim(tar.lm)[1])
                {n <-  dim(tar.lm)[1]
               }
            if (n == 0)
              {surf <- 1:lmdim
               ntmp <- lmdim
             }
            else
                { ntmp <- n
                  surf <- c(1:lmdim,sample((lmdim+1):dim(tar.lm)[1])[1:(n-lmdim)])
                }
            tar.lm <- tar.lm[surf,]
            ref.lm <- ref.lm[surf,]
            tar.norm <-  tar.norm[surf,]
                                        # print(tar.norm)
            L <- CreateL(ref.lm)
            U <- calcTang_U_s(tar.lm,tar.norm,SMvector=1:ntmp,deselect=FALSE,surface=1:ntmp)
                slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
                pro <- proj.read(slide,mesh2)
                
                tar.lm <- t(pro$vb[1:3,])
          }
        
### end relaxation
        if (!uselm)
          {tar.lm <- tar.lm[-c(1:lmdim),]
           ref.lm <- ref.lm[-c(1:lmdim),]
         }
      print(c(dim(ref.lm),dim(tar.lm)))
  
       tmp.mesh <- rotmesh.onto(tmp.mesh,ref.lm,tar.lm,adnormals=F)$mesh ## ro
       tmp <- rotmesh.onto(tmp.mesh1,ref.lm,tar.lm) ## rotate mesh

        tmp.mesh1 <- tmp$mesh
        tmp.lm <- tmp$yrot[1:lmdim,]
        
      
### check distance
        rs <- tmptar$quality
        rsme <- mean(tmptar$quality)

         }
      
     rgl.close()
   # dist <-  mean(sqrt(diag(tcrossprod(tar.lm-ref.lm))))
    return(list(mesh=tmp.mesh1,rsme=rsme,lm=tmp.lm,tar.lm=tar.lm))
   
  }
        
