#' mirror landmarks or triangular mesh in place
#'
#' mirror landmarks or triangular mesh in place
#'
#' @param x k x 3 matrix or mesh3d
#' @param icpiter integer: number of iterations to match reflected configuration onto original one
#' @param subsample integer: use only a subset for icp matching
#'
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
#' @export
mirror <- function(x,icpiter=50,subsample=NULL) UseMethod("mirror")

#' @rdname mirror
#' @export
mirror.matrix <- function(x,icpiter=50,subsample=NULL) {

    m <- ncol(x)
    if (m == 2)
        x <- cbind(x,0)
    
    xc <- apply(x,2,scale,scale=FALSE)
    trans <- x[1,]-xc[1,]
    pca <- prcomp(xc,scale. = F)
    krdelta <- diag(3)
    mirmat <- matrix(0,3,3)
    a <- c(0,0,1)
    anorm <- sum(a^2)
    for (i in 1:3)
        for (j in 1:3)
            mirmat[i,j] <- 2*a[i]*a[j]/anorm

    mirmat <- krdelta-mirmat
    out <- pca$x%*%t(mirmat)
    xrot = getTrafoRotaxis(c(0,0,0),c(1,0,0),pi)
    zrot = getTrafoRotaxis(c(0,0,0),c(0,0,1),pi)
    pca2 <- prcomp(out)
    test <- diag(pca2$rotation[,]%*%diag(3))
    if (test[3] < 0 || test[1] < 0)##test if rotation around x-axis is needed to fix orientation
        out <-  applyTransform(out,xrot)
    
    if (icpiter > 0)
        out <- icpmat(out,pca$x,icpiter,subsample = subsample)
    out <- out%*%t(pca$rotation)
    out <- t(t(out)+trans)
    if (m == 2)
        out <- out[,1:2]
    return(out)
    
}

#' @rdname mirror
#' @export
mirror.mesh3d <- function(x,icpiter=50,subsample=NULL) {
    mesh <- x
    x <- vert2points(mesh)
    vb <- mirror(x,icpiter=icpiter,subsample=subsample)
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
#' @param scale logical: if TRUE, scaling is allowed
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
icpmat <- function(x,y,iterations,mindist=1e15,subsample=NULL,scale=FALSE) {
    m <- ncol(x)
    if (m == 2) {
        x <- cbind(x,0)
        y <- cbind(y,0)
    }
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
        xtmp[,1:m] <- rotonmat(xtmp[,1:m],xtmp[good,1:m],y[clost$index[good],1:m],scale = scale,reflection=FALSE)
    }
    if (!is.null(subsample))
        xtmp <- rotonmat(x[,1:m],x[!subs,1:m],xtmp[,1:m],scale=scale,reflection=FALSE)
    if (m == 2)
        xtmp <- xtmp[,1:2]
    return(xtmp)
        
}
    
