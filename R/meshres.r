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
#' @export meshres
meshres <- function(mesh)
  {
      if (!inherits(mesh,"mesh3d"))
          stop("please provide object of class mesh3d")
      res <- 0
      VB <- mesh$vb[1:3,]
      nvb <- dim(VB)[2]
      IT <- mesh$it[1:3,]
      nit <- dim(IT)[2]
      storage.mode(VB) <- "double"
      storage.mode(IT) <- "integer"
      res <- .Fortran("meshres",VB,nvb,IT,nit,res)[[5]]
      return(res)
  }
    
