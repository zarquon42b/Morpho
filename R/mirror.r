#' mirror landmarks or triangular mesh in place
#'
#' mirror landmarks or triangular mesh in place
#'
#' @param x k x 3 matrix or mesh3d
#' @param icpiter integer: number of iterations to match reflected configuration onto original one
#' @param subsample integer: use only a subset for icp matching
#' @param pcAlign if TRUE, the icp will be preceeded by an alignment of the principal axis (only used if icpiter > 0), currently only works for 3D data.
#' @param mc.cores use parallel processing to find best alignment to original shape.
#' @param mirroraxis integer: which axis to mirror at
#' @param initPC logical: if TRUE the data will be prealigned by its principal axes.
#' @param initCenter logical: if TRUE and \code{initPC=FALSE}, \code{x} will be translated to its centroid before mirroring.
#' @details reflect a mesh configuration at the plane spanned by its first 2 principal axis, then try to rigidily register the reflected configuration onto the original one using iterative closest point search to establish correspondences.
#' @return returns the reflected object
#' @examples
#' data(boneData)
#' boneMir <- mirror(boneLM[,,1],icpiter=50,mc.cores=2,mirroraxis=3)
#' ## 2D Example:
#' if (require(shapes)) {
#' gorfMir <- mirror(gorf.dat[,,1],mirroraxis=2,pcAlign=TRUE,icpiter = 0)
#' plot(gorfMir,asp = 1)
#' points(gorf.dat[,,1],col=3)
#' }
#' \dontrun{
#' ## now mirror a complete mesh
#' require(rgl)
#' skullMir <- mirror(skull_0144_ch_fe.mesh,icpiter=10,subsample = 30,
#'                    mc.cores=2,mirroraxis=3,pcAlign=TRUE)
#' ###compare result to original
#' wire3d(skull_0144_ch_fe.mesh,col=3)
#' wire3d(skullMir,col=2)
#' }
#' @rdname mirror
#' @importFrom rgl rotationMatrix
#' @export
mirror <- function(x,icpiter=50,subsample=NULL,pcAlign=FALSE, mirroraxis=1,initPC=TRUE,initCenter=TRUE, mc.cores=2) UseMethod("mirror")

#' @rdname mirror
#' @export
mirror.matrix <- function(x,icpiter=50,subsample=NULL,pcAlign=FALSE, mirroraxis=1,initPC=TRUE,initCenter=TRUE,mc.cores=2) {

    m <- ncol(x)
    if (m == 2) {
        x <- cbind(x,0)
        pcAlign <- FALSE
    }
    if (initPC) {
        pca <- prcompfast(x,scale. = F)
        pca$rotation <- cbind(rbind(pca$rotation,0),0)
        pca$rotation[4,4] <- 1
    } else if (initCenter){
        xcenter <- scale(x,scale = F)
        rotmat <- computeTransform(xcenter,x)
        pca <- list(x=xcenter,rotation=rotmat,center=rep(0,3))
    } else {
        pca <- list(x=x,rotation=diag(m+1),center=rep(0,3))
    }

    ## i.e. a reflection along the z axis
    mirmat <- diag(c(1,1,1))
    mirmat[mirroraxis,mirroraxis] <- -1
    out <- pca$x%*%t(mirmat)
    
    if (pcAlign)
        out <- pcAlign(out,pca$x,iterations=icpiter,subsample = subsample,mc.cores = mc.cores)
    else if (icpiter > 0)
        out <- icpmat(out,pca$x,icpiter,subsample = subsample)
    
    
    out <- applyTransform(out,pca$rotation,inverse = !initPC)
    out <- t(t(out)+pca$center)
    if (m == 2)
        out <- out[,1:2]
    return(out)
    
}

#' @rdname mirror
#' @export
mirror.mesh3d <- function(x,icpiter=50,subsample=NULL,pcAlign=FALSE,mirroraxis=1,initPC=TRUE, initCenter=TRUE,mc.cores=2) {
    mesh <- x
    x <- vert2points(mesh)
    vb <- mirror(x,icpiter=icpiter,subsample=subsample,pcAlign=pcAlign,mirroraxis=mirroraxis,initPC=initPC,initCenter=initCenter, mc.cores=mc.cores)
    mesh$vb[1:3,] <- t(vb)
    mesh <- invertFaces(mesh)
    return(mesh)    
}

