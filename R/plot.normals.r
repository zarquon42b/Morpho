plot.normals <- function(x,long=1,lwd=1)
  {
    if ( ! "mesh3d" %in% class(x))
      {stop("please provide object of class mesh3d")
     }

    if (is.null(x$normals))
      {
        x <- adnormals(x)
      }

    n.mesh <- list()
    lvb <- dim(x$vb)[2]
    vb <- x$vb
    vb.norm <- x$vb+long*rbind(x$normals,0)
    vb <- cbind(vb,vb.norm)
    it <- rbind(1:lvb,1:lvb,(1:lvb)+lvb)
    n.mesh$vb <- vb
    n.mesh$it <- it
    class(n.mesh) <- c("mesh3d","shape3d")
    n.mesh$primitivetype <- "triangle"
    wire3d(n.mesh,color=3,lwd=lwd)
  }
    
      
