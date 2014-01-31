#' calculates the barycenters for all faces of a triangular mesh
#' 
#' calculates the barycenters for all faces of a triangular mesh
#' 
#' 
#' @param mesh triangular mesh of class 'mesh3d'
#' @return k x 3 matrix of barycenters for all \code{k} faces of input mesh.
#' @seealso \code{\link{closemeshKD}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' bary <- barycenter(shortnose.mesh)
#' \dontrun{
#' ##visualize mesh
#' wire3d(shortnose.mesh)
#' # visualize barycenters
#' points3d(bary, col=2)
#' ## now each triangle is equipped with a point in its barycenter
#' }
#' @export
barycenter <- function(mesh)
  {
    vb <- mesh$vb[1:3,]
    nvb <- dim(vb)[2]
    if (!is.matrix(vb) || !is.numeric(vb))
        stop("vertices must be a numeric matrix")
    if (!is.null(mesh$it))
        it <- mesh$it-1
    else
        stop("mesh has no triangular faces")
    nit <- dim(it)[2]
    out <- .Call("barycenter",vb,it)
    return(out)
  }
