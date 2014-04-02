#' separate a 3D-pointcloud by a hyperplane
#' separate a 3D-pointcloud by a hyperplane
#'
#' @param pointcloud numeric n x 3 matrix
#' @param v1 numeric vector of length=3 specifying a point on the separating plane
#' @param v2 numeric vector of length=3 specifying a point on the separating plane
#' @param v3 numeric vector of length=3 specifying a point on the separating plane
#' @param upper logical specify whether the points above or below the plane are to be reported as TRUE.
#' @return logical vector of length n. Reporting for each point if it is above or below the hyperplane
#' @details
#' As above and below are specified by the normal calculated from \eqn{(v2-v1) \times (v3-v1)}{(v2-v1) x (v3-v1)}, where \eqn{\times}{x} denotes the vector crossproduct. This means the normal points "upward" when viewed from the positon where v1, v2 and v3 are arranged counter-clockwise. Thus, which side is "up" depends on the ordering of v1, v2 and v3.
#' @examples
#' data(nose)
#' v1 <- shortnose.lm[1,]
#' v2 <- shortnose.lm[2,]
#' v3 <- shortnose.lm[3,]
#' pointcloud <- vert2points(shortnose.mesh)
#' upper <- cutSpace(pointcloud, v1, v2, v3)
#' \dontrun{
#' points3d(pointcloud[upper,])
#' }
#' @export
cutSpace <- function(pointcloud,v1, v2, v3,upper=TRUE) {
    e1 <- v2-v1
    e2 <- v3-v1
    e1 <- e1/sqrt(sum(e1^2))
    e2 <- e2/sqrt(sum(e2^2))

    normal <- crossp(e1,e2)
    normal <- normal/sqrt(sum(normal^2))
    e2a <- crossp(e1,normal)
    e2a <- e2a/sqrt(sum(e2a^2))
    Ep <- cbind(e1,e2a)
    pointcloud0 <- sweep(pointcloud,2,v1)
    orthopro <- t(Ep%*%t(Ep)%*%t(pointcloud0))
    diff <- pointcloud0-orthopro
    ins <- t(normal)%*%t(diff)
    #inside <- ang(diff,normal)
    if (upper)
        upside <- ins >= 0
    else
        upside <- ins <= 0
    return(upside)
}
#' cut a mesh by a hyperplane and remove parts above/below that plane
#'
#' cut a mesh by a hyperplane and remove parts above/below that plane
#' @param mesh triangular mesh of class "mesh3d"
#' @param v1 numeric vector of length=3 specifying a point on the separating plane
#' @param v2 numeric vector of length=3 specifying a point on the separating plane
#' @param v3 numeric vector of length=3 specifying a point on the separating plane
#' @param keep.upper logical specify whether the points above or below the plane are should be kept
#' @details see \code{\link{cutSpace}} for more details.
#' @return mesh with part above/below hyperplane removed
#' @export
cutMeshPlane <- function(mesh, v1, v2, v3, keep.upper=TRUE) {
    pointcloud <- vert2points(mesh)
    upper <- cutSpace(pointcloud, v1, v2, v3,upper=keep.upper)
    outmesh <- rmVertex(mesh,which(upper),keep = TRUE)
    return(outmesh)
}
                        
