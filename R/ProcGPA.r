#' Workhorse function for procSym, responsible for Procrustes registration
#' 
#' Workhorse function for procSym, responsible for Procrustes registration
#' 
#' 
#' @param dat.array Input k x m x n real array, where k is the number of
#' points, m is the number of dimensions, and n is the sample size.
#' @param tol numeric: Threshold for convergence during iterative
#' superimpositioning.
#' @param scale logical: indicating if scaling is requested
#' @param CSinit logical: if TRUE, all configurations are initially scaled to
#' Unit Centroid Size.
#' @param silent logical: suppress output of elapsed time.
#' @param weights numeric vector: assign per landmark weights.
#' @param centerweight logical: if TRUE, the landmark configuration is scaled
#' according to weights during the rotation process, instead of being scaled to
#' the Centroid size.
#' @param reflection logical: allow reflections.
#' @param pcAlign logical: if TRUE, the shapes are aligned by the principal axis of the first specimen, otherwise the orientation of the first specimen is used.
#' @return returns a list with
#' \item{rotated }{k x m x n array of the rotated configurations}
#' \item{mshape }{sample meanshape}
#' @author Stefan Schlager
#' @seealso \code{\link{procSym}, \link{rotonto}}
#' @references Goodall C. 1991. Procrustes methods in the statistical analysis
#' of shape. Journal of the Royal Statistical Society. Series B. Statistical
#' Methodology 53:285-239.
#' 
#' Dryden IL, Mardia KV.  1998. Statistical shape analysis. John Wiley and
#' sons, Chichester.
#' 
#' @examples
#' 
#' data(boneData)
#' proc <- ProcGPA(boneLM, CSinit=TRUE, silent=TRUE)
#' #now we landmarks 5 - 9 double the weight as  the others
#' weights <- c(rep(1,4),rep(2,5),1)
#' proc.wt <- ProcGPA(boneLM, CSinit=TRUE, weights=weights, silent=TRUE)
#' 
#' @export
ProcGPA <- function(dat.array,tol=1e-5,scale=TRUE,CSinit=FALSE,silent=FALSE,weights=NULL,centerweight=FALSE, reflection=TRUE,pcAlign=TRUE)
{
    if (!is.null(weights))
        weights <- weights/sum(weights)

    t0 <- Sys.time()
    x <- dat.array
    p1 <- 1e10
    p2 <- p1	
    n <- dim(dat.array)[3]
    k <- dim(dat.array)[1]
    m <- dim(dat.array)[2]
    x1 <- gdif(dat.array)
    
    arr.list <- list(0)	
###rotation step ####
    for ( i in 1:n)
        arr.list[[i]] <- list(x[,,i],1)
    
    if (CSinit) {
        arr.list <- lapply(arr.list, function(x){
                               x[[1]] <- scale(x[[1]], scale=FALSE);
                               x[[1]] <- x[[1]]/sqrt(sum(x[[1]]^2));
                               return(list(x[[1]],x[[2]]))
                           })
    } else { 
        arr.list <- lapply(arr.list,function(x){
                               x[[1]] <- scale(x[[1]], scale=FALSE);
                               return(list(x[[1]],x[[2]]))
                           })
    }
    mshape <- x[,,1]
    if (centerweight && !is.null(weights)) {
        mcent <- apply(mshape*weights,2,sum)           
        mshape <- scale(mshape,scale=F,center=mcent)
    }
### align mean by principal axes ###	
    if (pcAlign)
        mshape <- pcAlign(mshape)
    
    while (p1 > tol) {
### rotation of all configs on current consensus ###		
        arr.list <- lapply(arr.list,function(x){
                               x[[1]] <- rot.proc(x[[1]], x=mshape, scale=F,
                                                  weights=weights,
                                                  centerweight=centerweight,
                                                  reflection=reflection);
                               return(list(x[[1]],x[[2]]))
                           })
        
        for( i in 1:n)
            x[,,i] <- arr.list[[i]][[1]]
        
        x2 <- gdif(x)
        p1 <- abs(x1-x2)
        x1 <- x2
        
        
### scale/rotate step ###	
        if (scale) {      
            for ( i in 1:n)
                arr.list[[i]] <- list(x[,,i],1)
            
            while (p2 > tol) {
                for( i in 1:n)
                    if (!is.null(weights))
                        x[,,i] <- arr.list[[i]][[1]]*weights
                    else
                        x[,,i] <- arr.list[[i]][[1]]
                eigc <- scaleproc(x)
                for ( i in 1:n)	
                    arr.list[[i]][[2]] <- eigc[i]
                
                arr.list <- lapply(arr.list,function(x){
                                       x[[1]] <- x[[1]]*x[[2]];
                                       return(list(x[[1]],x[[2]]))
                                   })         
### rotation of all configs on current consensus ###		
                arr.list <- lapply(arr.list,function(x){
                                       x[[1]] <- rot.proc(x[[1]],x=mshape,scale=F,
                                                          weights=weights,
                                                          centerweight=centerweight,
                                                          reflection=reflection);
                                       return(list(x[[1]],x[[2]]))
                                   })
### scale step ####
                for( i in 1:n)
                    x[,,i] <- arr.list[[i]][[1]]
                
                x2 <- gdif(x)
                p2 <- abs(x1-x2)
                x1 <- x2
            }
        }
        mshape <- arrMean3(x)
        if (CSinit) {
            msize <- cSize(mshape)
            mshape <- mshape/msize
            if (scale)
                x <- x/msize
        }
    }
    t1 <- Sys.time()
    if (!silent)
        cat(paste("in... ",format(t1-t0)[[1]],"\n"))
    
    return(list(rotated=x,mshape=mshape))
}	
