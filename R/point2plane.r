#' projects a 3D coordinate orthogonally onto a plane
#'
#' projects a 3D coordinate orthogonally onto a plane
#' @param x 3D-vector or a k x 3 matrix with 3D vectors stored in rows
#' @param pt point on plane
#' @param pNorm plane normal (overrides specification by pt1 and pt2)
#' @param pt1 if pNorm=NULL, the plane will be defined by three points \code{pt, pt1, pt2}
#' @param pt2 if pNorm=NULL, the plane will be defined by three points \code{pt, pt1, pt2}
#' @return projected point
#' @examples
#' data(boneData)
#' ##project rhinion onto plane spanned by Nasion and both Nariales
#' rpro <- point2plane(boneLM[10,,1],pt=boneLM[9,,1],pt1=boneLM[3,,1],pt2=boneLM[4,,1])
#'
#' \dontrun{
#' require(rgl)
#' #visualize
#' wire3d(skull_0144_ch_fe.mesh,col="white")
#' ##get plane normal
#' normal <- crossProduct(boneLM[3,,1]-boneLM[9,,1],boneLM[4,,1]-boneLM[9,,1])
#' #' ## get plane offset
#' d <- norm(point2plane(c(0,0,0),pt=boneLM[9,,1],normal=normal),"2")
#' spheres3d(boneLM[,,1],radius=0.5)
#' spheres3d(boneLM[c(3,4,9),,1],radius=0.6,col=3)
#' ##original position of Rhinion
#' spheres3d(boneLM[10,,1],radius=0.6,col=2)
#' ##projected onto plane
#' spheres3d(rpro,radius=0.9,col=6)
#' lines3d(rbind(rpro,boneLM[10,,1]),lwd=3)
#' ##plot plane
#' planes3d(normal[1],normal[2],normal[3],d=d,col=2,alpha=0.5)
#'
#' ##now we project all points onto that plane:
#' spheres3d(point2plane(boneLM[,,1],pt=boneLM[9,,1],pt1=boneLM[3,,1],pt2=boneLM[4,,1]),col=3)
#'
#' ## and finally project the vertices of the mesh onto the plane
#' meshpro <- point2plane(vert2points(skull_0144_ch_fe.mesh),pt=boneLM[9,,1],normal=normal)
#' points3d(meshpro,col=2)
#' }
#' @rdname point2plane
#' @export
point2plane <- function(x, pt, normal=NULL, pt1=NULL, pt2=NULL) {
    if (is.vector(x))
        x <- matrix(x,1,3)

        if (is.null(normal) && is.null(pt1) && is.null(pt2))
        stop("either specify normal or pt1 and pt2")
    if (is.null(normal)) {
        e1 <- pt1-pt
        e1 <- e1/norm(e1,"2")
        e2 <- pt2-pt
        normal <- crossProduct(e1,e2)
        e2 <- crossProduct(e1,normal)
    } else {
        tp <- tangentPlane(normal)
        e1 <- tp$z
        e2 <- tp$y
    }
    Ep <- cbind(e1,e2)
    pointcloud0 <- sweep(x,2,pt)
    orthopro <- t(Ep%*%t(Ep)%*%t(pointcloud0))
    orthopro <- sweep(orthopro,2,-pt)
    if (nrow(orthopro) == 1)
        orthopro <- as.vector(orthopro)
    return(orthopro)
}
