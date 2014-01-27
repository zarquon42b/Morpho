#' calculate average edge length of a triangular mesh
#' 
#' calculate average edge length of a triangular mesh, by iterating over all
#' faces.
#' 
#' 
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @return returns average edge length (a.k.a. mesh resolution)
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
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
      if (!is.matrix(vb) || !is.numeric(vb))
          stop("vertices must be a numeric matrix")
      vb <- mesh$vb[1:3,]
      res <- .Call("meshres",vb,it)
      return(res)
  }
    
