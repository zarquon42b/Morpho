#' invert faces' orientation of triangular mesh
#' 
#' inverts faces' orientation of triangular mesh and recomputes vertex normals
#' 
#' 
#' @param mesh triangular mesh of class \code{mesh3d}
#' @return returns resulting mesh
#' @author Stefan Schlager
#' @seealso \code{\link{updateNormals}}
#' 
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' \dontrun{
#' shade3d(shortnose.mesh,col=3)
#' }
#' noseinvert <- invertFaces(shortnose.mesh)
#' ## show normals
#' \dontrun{
#' plotNormals(noseinvert,long=0.01)
#' }
#' @export
invertFaces <- function(mesh)
{ 	
	mesh$it <- mesh$it[c(3,2,1),]
        mesh <- vcgUpdateNormals(mesh)
  	return(mesh)
}
