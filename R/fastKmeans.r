#' fast kmeans clustering for 2D or 3D point clouds
#'
#' fast kmeans clustering for 2D or 3D point clouds - with the primary purpose to get a spatially equally distributed samples
#' @param x matrix containing coordinates or mesh3d
#' @param k number of clusters
#' @param iter.max maximum number of iterations
#' @param project logical: if x is a triangular mesh, the centers will be projected onto the surface.
#' @param threads integer number of threads to use
#' @return
#' returns a list containing
#' \item{selected}{coordinates closest to the final centers}
#' \item{centers}{cluster center}
#' \item{class}{vector with cluster association for each coordinate}
#' @examples
#' require(Rvcg)
#' data(humface)
#' set.seed(42)
#' clust <- fastKmeans(humface,k=1000,threads=1)
#' \dontrun{
#' require(rgl)
#'
#' ## plot the cluster centers
#' spheres3d(clust$centers)
#'
#' ## now look at the vertices closest to the centers
#' wire3d(humface)
#' spheres3d(vert2points(humface)[clust$selected,],col=2)
#' }
#' 
#' 
#' @export
fastKmeans <- function(x,k,iter.max=10,project=TRUE,threads=0) {
    isMesh <- FALSE
    if (inherits(x,"mesh3d")) {
        xorig <- x
        x <- vert2points(x)
        if (!is.null(xorig$it))
            isMesh <- TRUE
    }
    if (is.vector(x))
        x <- as.matrix(x)
    if (!ncol(x) %in% 1:3)
        stop("x can only have 1,2 or 3 columns")
    origdim <- ncol(x)
    if (origdim < 3) {
        supplement <- matrix(0,nrow(x),3-origdim)
        x <- cbind(x,supplement)
    }
    k <- abs(k)
    if (k >= nrow(x))
        stop("k exceeds sample size")

    centerinit <- sample(1:nrow(x))[1:k]
    centers <- x[centerinit,]
    cnt <- 1
    centerchk <- 1e12
    clost <- 1:nrow(x)
    while (cnt < iter.max && centerchk > 0) {
        indexold <- clost
        clost <- vcgKDtree(centers,x,k=1,threads = threads)$index
        centers <- .Call("fastSubsetMeans",x,clost-1L,k,threads)
        ## precaution for empty centers due to bad initialization
        if (sum(centers$checkempty)) { 
            centerinit <- sample(1:nrow(x))[1:k]
            centers <- x[centerinit,]
            cnt <- 1
        } else {
            centers <- centers$centers
        }
        centerchk <- max(abs(clost-indexold))
        
        cnt <- cnt+1
    }
    if (isMesh && project)
        centers <- vert2points(vcgClost(centers,xorig))
    clost_center <- sort(unique(vcgKDtree(x,centers,k=1,threads = threads)$index))
    out <- list(selected=clost_center,centers=centers[,1:origdim],class=clost)
    return(out)
}
    
