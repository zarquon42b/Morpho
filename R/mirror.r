#' mirror landmarks or triangular mesh in place
#'
#' mirror landmarks or triangular mesh in place
#'
#' @param x k x 3 matrix or mesh3d
#' @param icpiter integer: number of iterations to match reflected configuration onto original one
#' @param subsample integer: use only a subset for icp matching
#' @param pcAlign if TRUE, the icp will be preceeded by an alignment of the principal axis (only used if icpiter > 0).
#' @param mc.cores use parallel processing to find best alignment to original shape.
#' @details reflect a mesh configuration at the plane spanned by its first 2 principal axis, then try to rigidily register the reflected configuration onto the original one using iterative closest point search to establish correspondences.
#' @return returns the reflected object
#' @examples
#' data(boneData)
#' boneMir <- mirror(boneLM[,,1],icpiter=50)
#' ## 2D Example:
#' require(shapes)
#' gorfMir <- mirror(gorf.dat[,,1])
#' plot(gorfMir,asp = 1)
#' points(gorf.dat[,,1],col=3)
#' \dontrun{
#' ## now mirror a complete mesh
#' require(rgl)
#' skullMir <- mirror(skull_0144_ch_fe.mesh,icpiter=10,subsample = 30)
#' ###compare result to original
#' wire3d(skull_0144_ch_fe.mesh,col=3)
#' wire3d(skullMir,col=2)
#' }
#' @rdname mirror
#' @importFrom rgl rotationMatrix
#' @export
mirror <- function(x,icpiter=50,subsample=NULL,pcAlign=TRUE,mc.cores=2) UseMethod("mirror")

#' @rdname mirror
#' @export
mirror.matrix <- function(x,icpiter=50,subsample=NULL,pcAlign=TRUE,mc.cores=2) {

    m <- ncol(x)
    if (m == 2)
        x <- cbind(x,0)
    
    pca <- prcomp(x,scale. = F)
    # i.e. a reflection along the z axis
    mirmat <- diag(c(1,1,-1))
    out <- pca$x%*%t(mirmat)
    
    
    if (pcAlign)
        out <- pcAlign(out,pca$x,iterations=icpiter,subsample = subsample,mc.cores = mc.cores)
    else if (icpiter > 0)
        out <- icpmat(out,pca$x,icpiter,subsample = subsample)
    
    out <- out%*%t(pca$rotation)
    out <- t(t(out)+pca$center)
    if (m == 2)
        out <- out[,1:2]
    return(out)
    
}

#' @rdname mirror
#' @export
mirror.mesh3d <- function(x,icpiter=50,subsample=NULL,pcAlign=TRUE,mc.cores=2) {
    mesh <- x
    x <- vert2points(mesh)
    vb <- mirror(x,icpiter=icpiter,subsample=subsample,pcAlign=pcAlign,mc.cores=mc.cores)
    mesh$vb[1:3,] <- t(vb)
    mesh <- invertFaces(mesh)
    return(mesh)    
}

#' match two landmark configurations using iteratively closest point search
#'
#' match two landmark configurations using iteratively closest point search
#'
#' @param x moving landmarks
#' @param y target landmarks
#' @param iterations integer: number of iterations
#' @param mindist restrict valid points to be within this distance
#' @param subsample use a subsample determined by kmean clusters to speed up computation
#' @param type character: select the transform to be applied, can be "rigid","similarity" or "affine"
#' @return returns the rotated landmarks
#' @examples
#' data(nose)
#' icp <- icpmat(shortnose.lm,longnose.lm,iterations=10,subsample = 20)
#' 
#' ##2D example  using icpmat to determine point correspondences
#' require(shapes)
#' ## we scramble rows to show that this is independent of point order
#' moving <- gorf.dat[sample(1:8),,1]
#' plot(moving,asp=1) ## starting config
#' icpgorf <- icpmat(moving,gorf.dat[,,2],iterations = 20)
#' points(icpgorf,asp = 1,col=2)
#' points(gorf.dat[,,2],col=3)## target
#'
#' ## get correspondences using nearest neighbour search
#' index <- mcNNindex(icpgorf,gorf.dat[,,2],k=1,cores=1)
#' icpsort <- icpgorf[index,]
#' for (i in 1:8)
#' lines(rbind(icpsort[i,],gorf.dat[i,,2]))
#' @importFrom Rvcg vcgKDtree
#' @export
icpmat <- function(x,y,iterations,mindist=1e15,subsample=NULL,type=c("rigid","similarity","affine")) {
    m <- ncol(x)
    if (m == 2) {
        x <- cbind(x,0)
        y <- cbind(y,0)
    }
    type <- match.arg(type,c("rigid","similarity","affine"))
    if (!is.null(subsample)) {
        subsample <- min(nrow(x)-1,subsample)
        subs <- duplicated(kmeans(x,centers=subsample,iter.max =100)$cluster)
        xtmp <- x[!subs,]
    } else {
        xtmp <- x
    }
    for (i in 1:iterations) {
        clost <- vcgKDtree(y,xtmp,1)
        good <- which(clost$distance < mindist)
        trafo <- computeTransform(y[clost$index[good],],xtmp[good,],type=type)
        xtmp <- applyTransform(xtmp[,],trafo)
    }
    if (!is.null(subsample)) {
        fintrafo <- computeTransform(xtmp[,],x[!subs,],type = type)
        xtmp <- applyTransform(x,fintrafo)
    }
    if (m == 2)
        xtmp <- xtmp[,1:2]
    return(xtmp)
        
}
    
