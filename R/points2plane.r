#' projects a 3D coordinate orthogonally onto a plane
#'
#' projects a 3D coordinate orthogonally onto a plane
#' @param x 3D-vector or a k x 3 matrix with 3D vectors stored in rows
#' @param v1 point on plane
#' @param normal plane normal (overrides specification by v2 and v3)
#' @param v2 if pNorm=NULL, the plane will be defined by three points \code{v1, v2, v3}
#' @param v3 if pNorm=NULL, the plane will be defined by three points \code{v1, v2, v3}
#' @return projected point
#' @examples
#' data(boneData)
#' ##project rhinion onto plane spanned by Nasion and both Nariales
#' rpro <- points2plane(boneLM[10,,1],v1=boneLM[9,,1],v2=boneLM[3,,1],v3=boneLM[4,,1])
#'
#' \dontrun{
#' require(rgl)
#' #visualize
#' wire3d(skull_0144_ch_fe.mesh,col="white")
#' ##get plane normal
#' normal <- crossProduct(boneLM[3,,1]-boneLM[9,,1],boneLM[4,,1]-boneLM[9,,1])
#' #' ## get plane offset
#' d <- norm(points2plane(c(0,0,0),v1=boneLM[9,,1],normal=normal),"2")
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
#' spheres3d(points2plane(boneLM[,,1],v1=boneLM[9,,1],v2=boneLM[3,,1],v3=boneLM[4,,1]),col=3)
#'
#' ## and finally project the vertices of the mesh onto the plane
#' meshpro <- points2plane(vert2points(skull_0144_ch_fe.mesh),v1=boneLM[9,,1],normal=normal)
#' points3d(meshpro,col=2)
#' }
#' @rdname points2plane
#' @export
points2plane <- function(x, v1, normal=NULL, v2=NULL, v3=NULL) {
    if (is.vector(x))
        x <- matrix(x,1,3)

        if (is.null(normal) && is.null(v2) && is.null(v3))
        stop("either specify normal or v2 and v3")
    if (is.null(normal)) {
        e1 <- v2-v1
        e1 <- e1/norm(e1,"2")
        e2 <- v3-v1
        normal <- crossProduct(e1,e2)
        e2 <- crossProduct(e1,normal)
    } else {
        tp <- tangentPlane(normal)
        e1 <- tp$z
        e2 <- tp$y
    }
    Ep <- cbind(e1,e2)
    pointcloud0 <- sweep(x,2,v1)
    orthopro <- t(Ep%*%t(Ep)%*%t(pointcloud0))
    orthopro <- sweep(orthopro,2,-v1)
    if (nrow(orthopro) == 1)
        orthopro <- as.vector(orthopro)
    return(orthopro)
}

#' mirror points or mesh on an arbitrary plane
#'
#' mirror points or mesh on an arbitrary plane
#' @param x x 3D-vector or a k x 3 matrix with 3D vectors stored in rows. Or a triangular mesh of class mesh3d
#' @param v1 point on plane
#' @param normal plane normal (overrides specification by v2 and v3)
#' @param v2 if pNorm=NULL, the plane will be defined by three points \code{v1, v2, v3}
#' @param v3 if pNorm=NULL, the plane will be defined by three points \code{v1, v2, v3}
#' @return mirrored coordinates mesh
#' @examples
#' # mirror mesh on plane spanned by 3 midsagital landmarks
#' data(boneData)
#' mirrmesh <- mirror2plane(skull_0144_ch_fe.mesh,v1=boneLM[1,,1],v2=boneLM[9,,1],v3=boneLM[10,,1])
#' @rdname mirror2plane
#' @export
mirror2plane <- function(x,v1, normal=NULL, v2=NULL, v3=NULL) UseMethod("mirror2plane")

#'@rdname mirror2plane
#' @export
mirror2plane.matrix <- function(x,v1, normal=NULL, v2=NULL, v3=NULL){
    onplane <- points2plane(x,v1=v1,normal=normal,v2=v2,v3=v3)
    diff <- onplane-x
    mirrored <- onplane+diff
    return(mirrored)
}

#' @rdname mirror2plane
#' @export
mirror2plane.mesh3d <- function(x,v1, normal=NULL, v2=NULL, v3=NULL){
    mesh <- x
    x <- vert2points(x)
    newx <- mirror2plane(x,v1=v1,normal=normal,v2=v2,v3=v3)
    mesh$vb[1:3,] <- t(newx)
    mesh <- invertFaces(mesh)
    return(mesh)
}
