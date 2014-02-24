#' Estimate the shape by averaging the shape of the nearest neighbours.
#' 
#' Estimate the shape of one set of landmarks by averaging the shape of the
#' nearest neighbours obtained by a second set of landmarks. Weights are
#' calculated either form Mahalanobis or Procrustes distances. This can be
#' useful for data with missing landmarks.
#' 
#' This function calculates weights from one set of shape data and then
#' estimates the shape of another (or same) set of landmarks.  CAUTION:
#' landmark data has to be registered beforehand.
#' 
#' @param x an array or matrix (one row per specim) with data used for
#' estimating weights.
#' @param y an array or matrix (one row per specim) with landmark data on which
#' the weighted averaging is applied for prediction. If NULL, x will be used
#' for both tasks.
#' @param n amount of nearest neighbours to consider
#' @param mahalanobis logical: use mahalanobis distance
#' @param mc.cores integer: amount of cores used for parallel processing.
#' @return matrix or array of estimates.
#' @seealso \code{\link{proc.weight}}, \code{\link{fixLMtps}}
#' 
#' @examples
#' 
#' library(shapes)
#' proc <- procSym(gorf.dat)
#' #use the closest 3 specimen based on the first 4 landmarks
#' #to estimate the shape
#' estim <- NNshapeReg(proc$rotated[1:4,,],proc$rotated,n=3,mc.cores=1)
#' #compare estimation and true config
#' plot(proc$rotated[,,1],asp=1)
#' points(estim[,,1],col=2)
#' 
#' 
#' @export
NNshapeReg <- function(x,y=NULL, n=3, mahalanobis=FALSE,mc.cores = parallel::detectCores()) {
    if (is.null(y))
        y <- x
    outdim <- dim(y)
    if (length(dim(x)) == 3)
        x <- vecx(x)
    if (length(dim(y)) == 3)
        y <- vecx(y)
    i <- NULL
    win <- FALSE
    if(.Platform$OS.type == "windows")
        win <- TRUE
    else
        registerDoParallel(cores=mc.cores)### register parallel backend
    out <- y
    estfun <- function(i) {
        weighcalc <- proc.weight(x,n,i,mahalanobis=mahalanobis,report=F)$data
        ws <- diag(weighcalc$weight)
        tmpres <- apply(t(t(y[weighcalc$nr,])%*%ws),2,sum)
        return(tmpres)
    }
    if (win)
        out <- foreach(i=1:dim(x)[1],.combine=rbind) %do% estfun(i)
    else
        out <- foreach(i=1:dim(x)[1],.combine=rbind) %dopar% estfun(i)
    
    if (length(outdim) == 3) {
        out1 <- array(NA, dim=outdim)
        for (i in 1:outdim[3])
            out1[,,i] <- matrix(out[i,],outdim[1],outdim[2])
        
        out <- out1
    }
    return(out)
}

