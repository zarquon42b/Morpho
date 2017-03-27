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
#' @param weights vector of length \code{nrow(x)} containing weights for each row in \code{x}
#' @param centerweight logical: if weights are defined and centerweigths=TRUE, the matrix will be centered according to these weights instead of the barycenter.
#' @param threads integer: number of threads to use.
#' @return returns the rotated landmarks
#' @examples
#' data(nose)
#' icp <- icpmat(shortnose.lm,longnose.lm,iterations=10)
#'
#' ## example using weights
#' ## we want to assign high weights to the first three cordinates
#' icpw <- icpmat(shortnose.lm,longnose.lm,iterations=10,
#'                weights=c(rep(100,3),rep(1,620)),centerweight = TRUE)
#' ## the RMSE between those four points and the target is now smaller:
#' require(Rvcg)
#' RMSE <- sqrt(sum(vcgKDtree(longnose.lm,icp[1:3,],k=1)$distance^2))
#' RMSEW<- sqrt(sum(vcgKDtree(longnose.lm,icpw[1:3,],k=1)$distance^2))
#' barplot(c(RMSE,RMSEW),names.arg=c("RMSE weighted","RMSE unweighted"))
#' \dontrun{
#' ## plot the differences between unweighted and weighted icp
#' deformGrid3d(icp,icpw)
#' ## plot the first four coordinates from the icps:
#' spheres3d(icp[1:3,],col="red",radius = 0.5)
#' spheres3d(icpw[1:3,],col="green",radius = 0.5)
#' ## plot the target
#' spheres3d(longnose.lm,col="yellow",radius = 0.2)
#' }
#' ##2D example  using icpmat to determine point correspondences
#' if (require(shapes)) {
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
#' }
#' @importFrom Rvcg vcgKDtree vcgSearchKDtree vcgCreateKDtree
#' @export
icpmat <- function(x,y,iterations,mindist=1e15,subsample=NULL,type=c("rigid","similarity","affine"),weights=NULL,threads=1,centerweight=FALSE) {
    m <- ncol(x)
    if (m == 2) {
        x <- cbind(x,0)
        y <- cbind(y,0)
    }
    if (!is.null(weights))
        if (length(weights) != nrow(x))
            stop("weights must be of same length as nrow(x)")
    type <- match.arg(type,c("rigid","similarity","affine"))
    if (!is.null(subsample)) {
        subsample <- min(nrow(x)-1,subsample)
        subs <- fastKmeans(x,k=subsample,iter.max = 100,threads=threads)$selected
        xtmp <- x[subs,]
        if (!is.null(weights))
            weights <- weights[subs]
    } else {
        xtmp <- x
    }
    yKD <- vcgCreateKDtree(y)
    for (i in 1:iterations) {
        clost <- vcgSearchKDtree(yKD,xtmp,1,threads=threads)
        good <- which(clost$distance < mindist)
        tmpweights <- weights
        if (!is.null(weights))
            tmpweights <- weights[good]
        trafo <- computeTransform(y[clost$index[good],],xtmp[good,],type=type,weights = tmpweights,centerweight = centerweight)
        xtmp <- applyTransform(xtmp[,],trafo)
    }
    if (!is.null(subsample)) {
        fintrafo <- computeTransform(xtmp[,],x[subs,],type = type)
        xtmp <- applyTransform(x,fintrafo)
    }
    if (m == 2)
        xtmp <- xtmp[,1:2]
    return(xtmp)
        
}
    







