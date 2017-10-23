#' plots the normals of a triangular surface mesh.
#' 
#' visualises the vertex normals of a triangular surface mesh of class mesh3d.
#' If no normals are contained, they are computed.
#' 
#' 
#' @param x object of class "mesh3d"
#' @param length either a single numeric value or a numeric vector defining per-normals lenght (default is 1)
#' @param lwd width of the normals
#' @param col color of the normals
#' @param ... addtional parameters, currently not in use.
#' @author Stefan Schlager
#' 
#' @examples
#'
#' \dontrun{
#' require(rgl)
#' data(nose)
#' plotNormals(shortnose.mesh,col=4,long=0.01)
#' shade3d(shortnose.mesh,col=3)
#' }
#' 
#' @export
plotNormals <- function(x,length=1,lwd=1,col=1,...) {
    if ( ! "mesh3d" %in% class(x))
        stop("please provide object of class mesh3d")

    args <- list(...)
    print(args)
    if("long" %in% names(args)) {
        length <- args$long
        warning("argument 'long' is deprecated, please use 'length' instead")
    }
    
    if (is.null(x$normals)) {
        if (!is.null(x$it))
            x <- vcgUpdateNormals(x)
        else
            stop("mesh has neither normals nor faces")
    }

    n.mesh <- list()
    lvb <- dim(x$vb)[2]
    vb <- x$vb
    vb.norm <- vb[1:3,,drop=FALSE]+t(length*t(x$normals[1:3,,drop=FALSE]))
    vb <- cbind(vb[1:3,,drop=FALSE],vb.norm)
    vb <- rbind(vb,1)
    
    it <- rbind(1:lvb,1:lvb,(1:lvb)+lvb)
    n.mesh$vb <- vb
    n.mesh$it <- it
    class(n.mesh) <- c("mesh3d","shape3d")
   # n.mesh$primitivetype <- "triangle"
    wire3d(n.mesh,color=col,lwd=lwd,lit=FALSE)
    
  }
    
      
