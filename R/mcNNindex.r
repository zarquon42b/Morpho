#' find nearest neighbours for point clouds
#' 
#' find nearest neighbours for point clouds by using algorithms from the ANN
#' library. This is just a wrapper for the function ann from the package
#' yaImpute, enabling parallel processing.
#' 
#' wraps the function \code{ann} from package 'yaImpute' to allow multicore
#' processing
#' 
#' @param target \code{k x m} matrix containing data which to search.
#' @param query \code{l x m} matrix containing data for which to search.
#' @param cores integer: amount of CPU-cores to be used. Speed benefits are
#' only relevant for \code{k > 20}
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
#' 
#' @export
mcNNindex <- function(target,query,cores=parallel::detectCores(),k=k,...)
    {
        if(.Platform$OS.type == "windows")
            cores <- 1
        if (nrow (query) < 1000)
            cores <- 1
        out <- NULL
        mclist <- list()
        nx <- dim(query)[1]
        iter <- floor(nx/cores)
        if (cores > 1) {
            for (i in 1:(cores-1))
                mclist[[i]] <- query[(1:iter)+((i-1)*iter),]
            
            mclist[[cores]] <- query[-c(1:((cores-1)*iter)),]
        }  else
            mclist[[1]] <- query
        tmpfun <- function(x,...)
            {
                ##tmp0 <- nn2(target,x,k=k,searchtype="priority",...)$nn.idx
                tmp0 <- ann(ref=target, target=x, k=k, search.type="priority",verbose=FALSE)$knnIndexDist[,1:k] ## ann function from package yaImpute
                return(tmp0)
            }
        tmp <- mclapply(mclist,tmpfun,mc.cores=cores)
        if (k > 1) {
            for (i in 1:cores)
                out <- rbind(out,tmp[[i]])
        } else {
            out <- unlist(tmp)
        }
        return(out)
    }
