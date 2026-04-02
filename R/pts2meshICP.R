#' Align a landmark configuration to a mesh using ICP 
#'
#' Align a landmark configuration/pointcloud to a mesh using Iterative Closest Point search
#'
#' 
#' @param x moving points
#' @param y target mesh
#' @param maxiter integer: maximum number of iterations
#' @param minchange double: smallest maximum change between ICP interations
#' @param mindist restrict valid points to be within this distance
#' @param subsample use only a subsample of points to speed up computation, if subsampling is desired specify parameter k for k-means clustering (number of subsampled points)
#' @param weights vector of length \code{nrow(x)} containing weights for each row in \code{x}
#' @param centerweight logical: if weights are defined and centerweigths=TRUE, the matrix will be centered according to these weights instead of the barycenter.
#' @param type character: select the transform to be applied, can be "rigid","similarity" or "affine"
#' @param threads integer: number of threads to use.
#' @param minscale double: for types "similarity" and "affine", minimum scale factor per iteration, prevents excessive shrinking in first iteration.
#' @param mincumscale double: for types "similarity" and "affine", minimum cumulative scale factor, prevents shrinking below threshold (e.g. 0.7), values above 1 will enforce size increase, should only be used if both landmark configurations are in the same scale.
#' @return landmark array: moving points after the alignment
#' @examples
#' data(nose)
#' ## Example using rigid transformations
#' shortnose.lm[,3] <- shortnose.lm[,3] + 10
#' shortnose.lm = rotaxis3d(shortnose.lm, pt1=c(0,1,0), theta=pi/5)
#' pts2meshICP(shortnose.lm, shortnose.mesh, type = "rigid", subsample=500)
#' ## Example using affine transformations
#' pts2meshICP(longnose.lm, shortnose.mesh, type = "affine")
#' @importFrom Rvcg vcgClostKD
#' @export
pts2meshICP <- function(x,
                   y,
                   maxiter=100,
                   mindist=1e15,
                   minchange=0.0001,
                   subsample=NULL,
                   type=c("rigid","similarity","affine"),
                   threads=1,
                   weights=NULL,
                   centerweight=FALSE,
                   minscale=0.9,
                   mincumscale=NULL) {
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
  cumscalef <- 1
  i=1
  while(i <= maxiter) {
    yclost <- vert2points( vcgClostKD(xtmp, y) )
    distances <- sqrt(rowSums((yclost - xtmp)^2))
    good <- which(distances < mindist)
    tmpweights <- weights
    if (!is.null(weights))
      tmpweights <- weights[good]
    trafo <- computeTransform(yclost[good,],xtmp[good,],type=type,weights = tmpweights,centerweight = centerweight)
    if(type!="rigid"){
      # calculate scale factor from transformation matrix
      scalef <- sqrt(sum(trafo[1:3,1]^2))
      if (scalef<minscale){
        trafo <- computeTransform(yclost[good,],xtmp[good,], type = "rigid",weights = tmpweights,centerweight = centerweight)
        scalef <- 1
      }
      cumscalef <- cumscalef * scalef
      if(!is.null(mincumscale)){
        if(cumscalef<mincumscale){
          type <- "rigid"
          warning("Cumulative scale factor below mincumscale, switching to rigid type")
          #scale up by inverse of cumscalef
          trafo <- matrix(0, 4, 4) 
          trafo[1,1] <- 1/cumscalef
          trafo[2,2] <- 1/cumscalef
          trafo[3,3] <- 1/cumscalef
          cumscalef <- 1 # reset cumscalef
        }
      }
    }
    xtmpold <- xtmp
    xtmp <- applyTransform(xtmp,trafo)
    eucldist <- sqrt(rowSums((xtmpold - xtmp)^2))
    if(max(eucldist) < minchange){
      i=maxiter
    }
    i = i + 1
  }
  if (!is.null(subsample)) {
    fintrafo <- computeTransform(xtmp[,],x[subs,],type = type)
    xtmp <- applyTransform(x,fintrafo)
  }
  if (m == 2)
    xtmp <- xtmp[,1:2]
  return(xtmp)
}