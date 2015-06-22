#' align two 3D-pointclouds/meshes by their principal axes
#'
#' align two 3D-pointclouds/meshes by their principal axes
#' @param x matrix or mesh3d
#' @param y matrix or mesh3d, if missing, x will be centered by its centroid and aligned by its princial axis.
#' @param optim logical if TRUE, the RMSE between reference and target will be minimized testing all possible axes alignments and (if iterations > 0) followed by a rigid ICP procedure.
#' @param subsample integer: use subsampled points to decrease computation time of optimization.
#' @param iterations integer: number of iterations for optimization (the higher the more accurate but also more time consuming).
#' @param mc.cores use parallel processing to find best alignment to original shape.
#' @return rotated and translated version of x to the center and principal axes of y.
#' @details \code{x} and \code{y} will first be centered and aligned by their PC-axes. If \code{optim=TRUE},all possible 8 ordinations of PC-axes will be tested and the one with the smallest RMSE between the transformed version of \code{x} and the closest points on \code{y} will be used. Then the rotated version of \code{x} is translated to the original center of mass of \code{y}.
#' @examples
#' data(boneData)
#' blm1 <- pcAlign(boneLM[,,1],boneLM[,,2])
#' \dontrun{
#' require(rgl)
#' spheres3d(boneLM[,,1])#original position
#' spheres3d(blm1,col=2)#aligned configuration
#' spheres3d(boneLM[,,2],col=3)#target
#' }
#' @rdname pcAlign
#' @importFrom Rvcg vcgKDtree
#' @export
pcAlign <- function(x,y,optim=TRUE,subsample=NULL,iterations=10,mc.cores=2) UseMethod("pcAlign")

#' @rdname pcAlign
#' @export
pcAlign.matrix <- function(x, y,optim=TRUE,subsample=NULL,iterations=10, mc.cores=2) {
    if (!missing(y)) {
        if (inherits(y,"mesh3d"))
            y <- vert2points(y)
        if (.Platform$OS.type == "windows")
            mc.cores <- 1
        pca1 <- prcomp(x, retx=FALSE)
        pca2 <- prcomp(y, retx=FALSE)
        x <- scale(x, scale=F)    
        y <- scale(y, scale=F)
        rotx <- pca1$rotation
        roty <- pca2$rotation
        chk <- det(rotx%*%roty)
        if (chk  < 0)
            rotx[,3] <- -rotx[,3]

        x <- x%*%rotx
        y <- y%*%roty
        rotlist <- list(
            arot=getTrafoRotaxis(pt1=c(1,0,0),pt2=c(0,0,0),theta=pi),
            brot=getTrafoRotaxis(pt1=c(0,1,0),pt2=c(0,0,0),theta=pi),
            crot=getTrafoRotaxis(pt1=c(0,0,1),pt2=c(0,0,0),theta=pi))
        tests <- as.matrix(expand.grid(c(1,0),c(1,0),c(1,0)))
        tmpfun <- function(x,rotlist){
            for (i in 1:3) {
                if (x[i] == 0)
                    rotlist[[i]] <- diag(4)
            }
            return(rotlist)
        }
        subs <- rep(FALSE,nrow(x))
        if (!is.null(subsample)) {
            subsample <- min(nrow(x)-1,subsample)
            subs <- duplicated(kmeans(x,centers=subsample,iter.max =100)$cluster)
        }
        
        dists <- 1e10
        fintrafo <- diag(4)
       
        optifun <- function(i) {
            #for (i in 1:8) {
                rottmp <- tmpfun(tests[i,],rotlist)
                trafotmp <- rottmp[[1]]%*%rottmp[[2]]%*%rottmp[[3]]
                xtmp <- applyTransform(x,trafotmp)
                if (iterations > 0) {
                    xtmp1 <- icpmat(xtmp,y,iterations=iterations,subsample=subsample)
                    trafoicp <- computeTransform(xtmp1,xtmp)
                    trafotmp <- trafoicp%*%trafotmp
                    disttmp <- mean(vcgKDtree(y,xtmp1[!subs,],k=1)$dist^2)
                } else {
                    disttmp <- mean(vcgKDtree(y,xtmp[!subs,],k=1)$dist^2)
                }
                out <- list(dist=disttmp,trafo=trafotmp)
            }
        
        if (optim) {
            checkit <- mclapply(1:8,optifun,mc.cores=mc.cores)
            testdists <- sapply(checkit,function(x) x <- x$dist)
            mindist <- which(testdists == min(testdists))[1]
            fintrafo <- checkit[[mindist]]$trafo
        }
        x <- applyTransform(x,fintrafo)
        x <- x%*%t(pca2$rotation)
        x <- t(t(x)+pca2$center)
                
        return(x)
    } else {
        x <- scale(x,scale=F)
        rotms <- eigen(crossprod(x))$vectors
        if (det(rotms) < 0)
            rotms[,1] <- rotms[,1]*-1
        
        x <- x%*%rotms
        return(x)
    }
}
#' @rdname pcAlign
#' @export
pcAlign.mesh3d <- function(x,y,optim=TRUE,subsample=NULL,iterations=10,mc.cores=2) {
    xorig <- x
    x <- vert2points(x)
    tmpverts <- pcAlign(x,y,optim=optim,subsample=subsample,iterations=iterations,mc.cores=mc.cores)
    xorig$vb[1:3,] <- t(tmpverts)
    xorig <- vcgUpdateNormals(xorig)
    return(xorig)
}
