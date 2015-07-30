#' separate a 3D-pointcloud by a hyperplane
#' 
#' separate a 3D-pointcloud by a hyperplane
#'
#' @param pointcloud numeric n x 3 matrix
#' @param v1 numeric vector of length=3 specifying a point on the separating plane
#' @param v2 numeric vector of length=3 specifying a point on the separating plane
#' @param v3 numeric vector of length=3 specifying a point on the separating plane
#' @param normal plane normal (overrides specification by v2 and v3)
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
#' require(rgl)
#' normal <- crossProduct(v2-v1,v3-v1)
#' zeroPro <- points2plane(rep(0,3),v1,normal)
#' ## get sign of normal displacement from zero
#' sig <- sign(crossprod(-zeroPro,normal))
#' d <- sig*norm(zeroPro,"2")
#' planes3d(normal[1],normal[2],normal[3],d=d)
#' points3d(pointcloud[upper,])
#' }
#' @export
cutSpace <- function(pointcloud,v1, v2=NULL, v3=NULL,normal=NULL, upper=TRUE) {
    orthopro <- points2plane(pointcloud,v1=v1,v2=v2,v3=v3,normal=normal)
    diff <- pointcloud-orthopro
    if (is.null(normal)) {
        e1 <- v2-v1
        e1 <- e1/norm(e1,"2")
        e2 <- v3-v1
        e2 <- e2/norm(e2,"2")
        normal <- crossProduct(e1,e2)
    }
    ins <- t(normal)%*%t(diff)
    #inside <- ang(diff,normal)
    if (upper)
        upside <- ins > 0
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
#' @param normal plane normal (overrides specification by v2 and v3)
#' @param keep.upper logical specify whether the points above or below the plane are should be kept
#' @details see \code{\link{cutSpace}} for more details.
#' @return mesh with part above/below hyperplane removed
#' @export
cutMeshPlane <- function(mesh, v1, v2=NULL, v3=NULL, normal=NULL,keep.upper=TRUE) {
    pointcloud <- vert2points(mesh)
    upper <- cutSpace(pointcloud, v1=v1,v2=v2,v3=v3,normal=normal,upper=keep.upper)
    outmesh <- list()
    lremain <- length(which(upper))
    if (lremain) {
        if (lremain < ncol(mesh$vb))
            outmesh <- rmVertex(mesh,which(upper),keep = TRUE)
        else
            outmesh <- mesh
    } else {
        warning("nothing left")
    }
    return(outmesh)
}
                        
