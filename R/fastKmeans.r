#' fast kmeans clustering for 2D or 3D point clouds
#'
#' fast kmeans clustering for 2D or 3D point clouds - with the primary purpose to get a spatially equally distributed samples
#' @param x matrix containing coordinates or mesh3d
#' @param k number of clusters
#' @param iter.max maximum number of iterations
#' @param tol convergence threshold if the mean distance between centers is below it, the algorithm stops
#' @param threads integer number of threads to use
#' @return
#' returns a list containing
#' \item{selected}{coordinates closest to the final centers}
#' \item{centers}{cluster center}
#' \item{class}{vector with cluster association for each coordinate}
#' @examples
#' require(Rvcg);require(rgl)
#' data(humface)
#' clust <- fastKmeans(humface,k=1000,threads=2)
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
fastKmeans <- function(x,k,iter.max=10,tol=1e-5,threads=parallel::detectCores()) {
    if (inherits(x,"mesh3d"))
        x <- vert2points(x)
    k <- abs(k)
    if (k >= nrow(x))
        return(list(centers=x,selected=1:nrow(x)))
    if (!ncol(x) %in% 2:3)
        stop("x can only have 2 or 3 columns")
    centerinit <- sample(1:nrow(x))[1:k]
    centers <- x[centerinit,]
    cnt <- 1
    centerchk <- 1e12
    while (cnt < iter.max && centerchk > tol) {
        centerold <- centers
        clost <- vcgKDtree(centers,x,k=1,threads = threads)$index
        centers <- .Call("fastSubsetMeans",x,clost-1L)
        centerchk <- mean(vcgKDtree(centers,centerold,k=1,threads = threads)$distance)
        clost_center <- sort(unique(vcgKDtree(x,centers,k=1,threads = threads)$index))
        cnt <- cnt+1
    }
    out <- list(selected=clost_center,centers=centers,class=clost)
    return(out)
}
    
