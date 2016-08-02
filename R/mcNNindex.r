#' find nearest neighbours for 2D and 3D point clouds
#' 
#' find nearest neighbours for point clouds using a kd-tree search. This is just a wrapper of the function vcgKDtree from
#' package Rvcg. Wwraps the function \code{vcgKDtree} from package 'Rvcg' (for backward compatibility )
#' 
#' @param target \code{k x m} matrix containing data which to search.
#' @param query \code{l x m} matrix containing data for which to search.
#' @param cores integer: amount of CPU-cores to be used. Only available on systems with OpenMP support.
#' @param k integer: how many closest points are sought.
#' @param \dots additional arguments - currently unused.
#' 
#' @return \code{l x k } matrix containing indices of closest points.
#' @seealso \code{\link{closemeshKD}}
#' 
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' # find closest vertex on surface for each landmark
#' clost <- mcNNindex(vert2points(shortnose.mesh),shortnose.lm, k=1,
#' mc.cores=1)
#' \dontrun{
#' spheres3d(vert2points(shortnose.mesh)[clost,],col=2,radius=0.3)
#' spheres3d(shortnose.lm,radius=0.3)
#' wire3d(shortnose.mesh)
#' }
#' @importFrom Rvcg vcgKDtree
#' @export
mcNNindex <- function(target,query,cores=parallel::detectCores(),k=k,...)
    {
        ## if(.Platform$OS.type == "windows")
        ##    cores <- 1
        ## if (nrow (query) < 1000)
        ##    cores <- 1
        ## out <- NULL
        ## mclist <- list()
        ## nx <- dim(query)[1]
        ## iter <- floor(nx/cores)
        ## if (cores > 1) {
        ##     for (i in 1:(cores-1))
        ##         mclist[[i]] <- query[(1:iter)+((i-1)*iter),]
            
        ##     mclist[[cores]] <- query[-c(1:((cores-1)*iter)),]
        ## }  else
        ##     mclist[[1]] <- query
        ## tmpfun <- function(x,...)
        ##     {
        ##         ##tmp0 <- nn2(target,x,k=k,searchtype="priority",...)$nn.idx
        ##         tmp0 <- ann(ref=target, target=x, k=k, search.type="priority",verbose=FALSE)$knnIndexDist[,1:k] ## ann function from package yaImpute
        ##         return(tmp0)
        ##     }
        ## tmp <- mclapply(mclist,tmpfun,mc.cores=cores)
        out <- vcgKDtree(target,query,k,threads=cores)$index
        
        return(out)
    }
