#' scale a mesh of class "mesh3d"
#' 
#' scales (the vertices of a mesh by a scalar
#' 
#' The mesh's center is determined either as mean of the bounding box
#' (center="bbox") or mean of vertex coordinates (center="mean") and then
#' scaled according to the scaling factor. If center="none", vertex coordinates
#' will simply be multiplied by "size".
#' 
#' @param mesh object of class "mesh3d"
#' @param size numeric: scale factor
#' @param center character: method to position center of mesh after scaling:
#' values are "bbox", and "mean". See Details for more info.
#' @return returns a scaled mesh
#' @author Stefan Schlager
#' @seealso \code{\link{rotmesh.onto}}
#' 
#' @examples
#' 
#' data(nose)
#' #inflate mesh by factor 4
#' largenose <- scalemesh(shortnose.mesh,4)
#' 
#' @export
scalemesh <- function(mesh,size,center=c("bbox","mean", "none"))
{	
    getmean <- TRUE
  if (substr(center[1],1L,1L) =="b")
    meshmean <- colMeans(meshcube(mesh))
  else if (substr(center[1],1L,1L) =="m")
    meshmean <- colMeans(vert2points(mesh))
  else if (substr(center[1],1L,1L) =="n")
     getmean <- FALSE
  else
    stop("Please provide valid centering method\n")
    if (getmean)
        {
            mesh <- translate3d(mesh,-meshmean[1],-meshmean[2],-meshmean[3])
            mesh$vb[1:3,] <- mesh$vb[1:3,]*size
            mesh <- translate3d(mesh,meshmean[1],meshmean[2],meshmean[3])
        }
   else
       mesh$vb[1:3,] <- mesh$vb[1:3,]*size
  return(mesh)
}
	
