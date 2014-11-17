#' projects a 3D coordinate orthogonally onto a plane
#'
#' projects a 3D coordinate orthogonally onto a plane
#' @param x 3D-vector
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
#' pNorm <- crossp(boneLM[3,,1]-boneLM[9,,1],boneLM[4,,1]-boneLM[9,,1])
#' pNorm <- pNorm/norm(pNorm,"2")
#' ## get plane offset
#' d <- norm(point2plane(c(0,0,0),pt=boneLM[9,,1],pt1=boneLM[3,,1],pt2=boneLM[4,,1]),"2")
#' spheres3d(boneLM[,,1],radius=0.5)
#' spheres3d(boneLM[c(3,4,9),,1],radius=0.6,col=3)
#' ##original position of Rhinion
#' spheres3d(boneLM[10,,1],radius=0.6,col=2)
#' ##projected onto plane
#' spheres3d(rpro,radius=0.9,col=6)
#' lines3d(rbind(rpro,boneLM[10,,1]),lwd=3)
#' ##plot plane
#' planes3d(pNorm[1],pNorm[2],pNorm[3],d=d,col=2,alpha=0.5)
#'
#' ##now we project all points onto that plane:
#' spheres3d(point2plane(boneLM[,,1],pt=boneLM[9,,1],pt1=boneLM[3,,1],pt2=boneLM[4,,1]),col=3)
#' }
#' @rdname point2plane
#' @export
#'
point2plane <- function(x, pt, pNorm=NULL, pt1=NULL, pt2=NULL) UseMethod("point2plane")

#' @rdname point2plane
#' @export
point2plane.default <- function(x, pt, pNorm=NULL, pt1=NULL, pt2=NULL) {
    if (is.null(pNorm) && is.null(pt1) && is.null(pt2))
        stop("either specify pNorm or pt1 and pt2")
    if (is.null(pNorm)) {
        v1 <- pt1-pt
        v2 <- pt2-pt
        pNorm <- crossp(v1,v2)
    }
    pNorm <- pNorm/norm(pNorm,"2")
    xdiff <- x-pt
    xnorm <- norm(xdiff,"2")
    ang <- angle.calc(xdiff,pNorm)
    if (is.nan(ang))
        ang <- 0
    if (ang > pi/2) {
        pNorm <- -pNorm
        ang <- pi-ang
    }
    nlen <- cos(ang)*xnorm
    xpro <- x-pNorm*nlen
    return(xpro)
}

#' @rdname point2plane
#' @export
point2plane.matrix <- function(x, pt, pNorm=NULL, pt1=NULL, pt2=NULL) {
    if (is.null(pNorm) && is.null(pt1) && is.null(pt2))
        stop("either specify pNorm or pt1 and pt2")
    if (is.null(pNorm)) {
        v1 <- pt1-pt
        v2 <- pt2-pt
        pNorm <- crossp(v1,v2)
    }
    pNorm <- pNorm/norm(pNorm,"2")
    xdiff <- sweep(x,2,pt)
    xnorm <- sqrt(rowSums(xdiff^2))
    ang <- apply(xdiff,1,angle.calc,y=pNorm)
    if (length(is.nan(ang)))
        ang[which(is.nan(ang))] <- 0
    sapply(ang,function(ang) {
        if (ang > pi/2) {
            pNorm <- -pNorm
            ang <- pi-ang
        } else
            ang <- ang
        return(ang)
    })
    nlen <- cos(ang)*xnorm
    xpro <- x-nlen*matrix(pNorm,nrow(x),ncol(x),byrow = T)
    return(xpro)
}
