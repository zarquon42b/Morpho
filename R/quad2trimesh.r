#' converts a mesh containing quadrangular faces into one only consisting of triangles
#'
#' converts a mesh containing quadrangular faces into one only consisting of triangles
#' @param mesh object of class "mesh3d"
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

quad2trimesh <- function(mesh) {
    if (!inherits(mesh,"mesh3d"))
        stop("please provide mesh of class mesh3d")
    if (is.null(mesh$ib)) {
        warning("this is no quadmesh, nothing to be done")
    } else {
        ib2it <- rbind(mesh$ib[1:3,],mesh$ib[c(3:4,1),])
        mesh$it <- ib2it
        mesh$ib <- NULL
        mesh <- updateNormals(mesh)
    }
    return(mesh)
}
