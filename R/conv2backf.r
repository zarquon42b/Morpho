#' invert faces' orientation of triangular mesh
#' 
#' inverts faces' orientation of triangular mesh and recomputes vertex normals
#' 
#' 
#' @param mesh triangular mesh of class \code{mesh3d}
#' @return returns resulting mesh
#' @author Stefan Schlager
#' @seealso \code{\link{adnormals}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' shade3d(shortnose.mesh,col=3)
#' noseinvert <- conv2backf(shortnose.mesh)
#' ## show normals
#' plotNormals(noseinvert,long=0.01)
#' 
#' @export conv2backf
conv2backf <- function(mesh)
{ 	
	mesh$it <- mesh$it[c(3,2,1),]
        mesh <- adnormals(mesh)
  	return(mesh)
}
