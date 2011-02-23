map.mesh <- function(mesh1,lm1,mesh2,lm2,tol=1e-3,it=2)
  {
    round <- 0
    p <- 1e10
    rsme <- 1e11
    tmp.mesh <- mesh1
    rot <- rotmesh.onto(mesh1,lm1,lm2)
    tmp.mesh <- rot$mesh
    ref.lm <- t(tmp.mesh$vb[1:3,])
    #tar.lm <- ref.lm <- t(tmp.mesh$vb[1:3,])
    
    while(p > tol && round < it )
      {
        round <- round+1
        rsme_old <- rsme
        tmp.mesh_old <- tmp.mesh
        tmptar <- proj.read(t(tmp.mesh$vb[1:3,]),mesh2)
        #quant <- quantile(tmptar$quality)
        #good <- which(tmptar$quality < quant[2])
        tar.lm <- t(tmptar$vb[1:3,])
        ref.lm <- t(tmp.mesh$vb[1:3,])
        tmp.mesh <- rotmesh.onto(tmp.mesh,ref.lm,tar.lm)$mesh
        #ref.lm <- t(tmp.mesh$vbood])

### check distance
        rs <- tmptar$quality
        rsme <- mean(tmptar$quality)

         }
      }
     gc()
   # dist <-  mean(sqrt(diag(tcrossprod(tar.lm-ref.lm))))
    return(list(mesh=tmp.mesh,rsme=rsme))
   
  }
        
