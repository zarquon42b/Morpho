
map.mesh <- function(mesh1,lm1,mesh2,lm2,tol=1e-3,it=2,overlap=0.8,raytol=NULL,strict=TRUE,n=NULL,uselm=TRUE,subset=NULL)
  {
    round <- 0
    p <- 1e10
    rsme <- 1e11
    adnormals <- TRUE
    cloud <- FALSE
    rot <- rotmesh.onto(mesh1,lm1,lm2)
    tmp.mesh1 <- rot$mesh

    if (!is.null(subset))
      if ("mesh3d" %in% class(subset))
        {
          rotsub <- rotmesh.onto(subset,lm1,lm2)
          tmp.mesh <- rotsub$mesh
        }
      else
### interactive selection of mesh region
        {
          adnormals <- FALSE
          cloud <- TRUE
          open3d()
          wire3d(mesh1,col=3)
          selcheck <- 0
          
          cat("select a region using the right mouse button\n")
          while (selcheck == 0)             
            {
             rgl.bringtotop(stay = FALSE)
              if (interactive())
                { f <-  select3d("right")
                  subset <- t(mesh1$vb[1:3,])
                  selected <- which(f(subset))
                  selcheck <- length(selected)
                  if (selcheck != 0)
                    {
                      subset <- subset[selected,]
                      view <- points3d(subset,col=2,cex=2)
                     
                      answer <- readline("do you like the view? (y/n)\n")
                      if (answer == "n")
                        {
                          selcheck <- 0
                          rgl.pop("shapes",id=view)   
                        }   
                    }
                  else
                    { cat("nothing selected")
                    }
                       
                }
               
              }
          ## end selection
          
          rotsub <- rotonmat(subset,lm1,lm2,scale=F)
          tmp.mesh <- proj.read(rotsub,tmp.mesh1)
          tmp.mesh$vb <- rbind(tmp.mesh$vb,1)
          rgl.close()
        }     
               
    tmp.lm <- rot$yrot
    ref.lm <- t(tmp.mesh$vb[1:3,])
    lmdim <- dim(lm1)[1]
    
    while(p > tol && round < it )
      {
        round <- round+1
        rsme_old <- rsme
        tmp.mesh_old <- tmp.mesh
        tmptar <- rbind(tmp.lm,t(tmp.mesh$vb[1:3,])) ### add landmarks to vertices of rotated reference mesh
        
        if (is.null(raytol))
          {
            ## project coordinates on target
            tmptar <- proj.read(tmptar,mesh2)
          }

        
        else

          {
            ## project coordinates on target along ray
            tmptar <- proj.read(tmptar,tmp.mesh1)           
            tmptar <- ray2mesh(tmptar,mesh2,tol=raytol,strict=strict)
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

### use a subset of vertices and relax against reference
        if (!is.null(n)) 
          {
            if (n > dim(tar.lm)[1])
              {
                n <-  dim(tar.lm)[1]
              }
            if (n == 0)
              {
                surf <- 1:lmdim
                ntmp <- lmdim
              }
            else
              {
                ntmp <- n
                surf <- c(1:lmdim,sample((lmdim+1):dim(tar.lm)[1])[1:(n-lmdim)])
              }
            tar.lm <- tar.lm[surf,]
            ref.lm <- ref.lm[surf,]
            tar.norm <-  tar.norm[surf,]
                                        
            L <- CreateL(ref.lm)
            U <- calcTang_U_s(tar.lm,tar.norm,SMvector=1:ntmp,deselect=FALSE,surface=1:ntmp)
            slide <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=3)$Gamatrix
            pro <- proj.read(slide,mesh2)
            
            tar.lm <- t(pro$vb[1:3,])
          }
### end relaxation
        
        if (!uselm)
          {
            tar.lm <- tar.lm[-c(1:lmdim),]
            ref.lm <- ref.lm[-c(1:lmdim),]
          }

        tmp <- rotmesh.onto(tmp.mesh1,ref.lm,tar.lm) ## rotate mesh
        tmp.mesh1 <- tmp$mesh
        tmp.lm <- tmp$yrot[1:lmdim,]

        if (cloud)
          {
            tmp.mesh <- rotmesh.onto(tmp.mesh,ref.lm,tar.lm,adnormals=adnormals)$mesh 
          }
        else
          {
            tmp.mesh <- tmp.mesh1
          }
        
### check distance
        rs <- tmptar$quality
        rsme <- mean(tmptar$quality)
        
      }
                                        # dist <-  mean(sqrt(diag(tcrossprod(tar.lm-ref.lm))))
    return(list(mesh=tmp.mesh1,rsme=rsme,lm=tmp.lm,tar.lm=tar.lm))
  }

