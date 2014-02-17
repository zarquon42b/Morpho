#' calculate average edge length of a triangular mesh
#' 
#' calculate average edge length of a triangular mesh, by iterating over all
#' faces.
#' 
#' 
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @return returns average edge length (a.k.a. mesh resolution)
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' data(boneData)
#' mres <- meshres(skull_0144_ch_fe.mesh)
#' 
#' 
#' @export
meshres <- function(mesh)
  {
      if (!inherits(mesh,"mesh3d"))
          stop("please provide object of class mesh3d")
      if (!is.null(mesh$it))
          it <- mesh$it-1
      else
          stop("mesh has no triangular faces")
      vb <- mesh$vb[1:3,]
       if (!is.matrix(vb) || !is.numeric(vb))
          stop("vertices must be a numeric matrix")
      res <- .Call("meshres",vb,it)
      return(res)
  }
    
