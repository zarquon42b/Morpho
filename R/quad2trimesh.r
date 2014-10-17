#' converts a mesh containing quadrangular faces into one only consisting of triangles
#'
#' converts a mesh containing quadrangular faces into one only consisting of triangles
#' @param mesh object of class "mesh3d"
#' @param updateNormals logical: request recalculation of (angle weighted) vertex normals.
#' @return triangular mesh with updated normals
#' @examples
#' 
#' Sigma <- diag(3:1) #create a 3D-covariance matrix
#' require(rgl)
#' quadmesh <- ellipse3d(Sigma)##create quadmesh
#' trimesh <- quad2trimesh(quadmesh)# convert to trimesh
#'
#'
#' @export

quad2trimesh <- function(mesh, updateNormals=TRUE) {
    if (!inherits(mesh,"mesh3d"))
        stop("please provide mesh of class mesh3d")
    if (is.null(mesh$ib)) {
        warning("this is no quadmesh, nothing to be done")
    } else {
        ib2it <- cbind(mesh$ib[1:3,,drop=FALSE],mesh$ib[c(3:4,1),,drop=FALSE])
        mesh$it <- ib2it
        mesh$ib <- NULL
        if (updateNormals) {
            mesh <- vcgUpdateNormals(mesh)
        }
    }
    return(mesh)
}
