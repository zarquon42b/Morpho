#' Project points onto the closest point on a mesh
#' 
#' project points onto a given surface and return projected points and normals.
#' 
#' 
#' @param lm m x 3 matrix containing 3D coordinates.
#' @param mesh character: specify path to mesh file.
#' @param readnormals logical: return normals of projected points.
#' @param smooth logical: rerturn smoothed normals.
#' @param sign logical: request signed distances.
#' @param \dots additional arguments currently not used.
#' @return if readnormals = FALSE, a m x 3 matrix containing projected points
#' is returned, otherwise a list, where
#' \item{vb }{3 x m matrix containing projected points}
#' \item{normals }{3 x m matrix containing normals} 
#' \item{quality }{vector containing distances }
#' @author Stefan Schlager
#' @seealso \code{\link{closemeshKD}}
#' @references Detection of inside/outside uses the algorithm proposed in:
#' 
#' Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance
#' Fields From Triangle Meshes. Informatics and Mathematical Modelling.
#' 
#' @examples
#' 
#' 
#' data(nose)
#' \dontrun{
#' repro <- projRead(shortnose.lm,shortnose.mesh)
#' }
#' 
#' @importFrom Rvcg vcgClost vcgImport
#' @export
projRead <- function(lm, mesh,readnormals=TRUE, smooth=FALSE, sign=TRUE,...)
{
    if (is.character(mesh))
        mesh <- vcgImport(mesh,updateNormals=FALSE,clean=FALSE)
    
    data <- vcgClost(lm, mesh, smoothNormals=smooth,sign=sign,borderchk=FALSE)
    if (!readnormals)
        data <- vert2points(data)
    return(data)
}
