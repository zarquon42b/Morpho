#' plots the normals of a triangular surface mesh.
#' 
#' visualises the vertex normals of a triangular surface mesh of class mesh3d.
#' If no normals are contained, they are computed.
#' 
#' 
#' @param x object of class "mesh3d"
#' @param long length of the normals (default is 1)
#' @param lwd width of the normals
#' @param col color of the normals
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' 
#' plotNormals(shortnose.mesh,col=4,long=0.01)
#' \dontrun{
#' shade3d(shortnose.mesh,col=3)
#' }
#' 
#' @export plotNormals
plotNormals <- function(x,long=1,lwd=1,col=1)
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
    vb.norm <- vb+long*rbind(x$normals[1:3,],0)
    vb.norm[4,] <- 1
    vb <- cbind(vb,vb.norm)
    it <- rbind(1:lvb,1:lvb,(1:lvb)+lvb)
    n.mesh$vb <- vb
    n.mesh$it <- it
    class(n.mesh) <- c("mesh3d","shape3d")
   # n.mesh$primitivetype <- "triangle"
    wire3d(n.mesh,color=col,lwd=lwd,lit=FALSE)
    
  }
    
      
