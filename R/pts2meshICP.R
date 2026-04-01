#' Align a landmark configuration to a mesh using ICP 
#'
#' Align a landmark configuration/pointcloud to a mesh using Iterative Closest Point search
#'
#' param weights vector of length \code{nrow(x)} containing weights for each row in \code{x}
#' param centerweight logical: if weights are defined and centerweigths=TRUE, the matrix will be centered according to these weights instead of the barycenter.
#' param threads integer: number of threads to use.
#' param subsample use only a subsample of points to speed up computation, if subsampling is desired specify parameter k for k-means clustering (number of subsampled points) 
#' 
#' @param x moving points
#' @param y target mesh
#' @param maxiter integer: maximum number of iterations
#' @param minchange double: smallest maximum change between ICP interations
#' @param mindist restrict valid points to be within this distance
#' @param type character: select the transform to be applied, can be "rigid","similarity" or "affine"
#' @param minscale double: for types "similarity" and "affine", minimum scale factor per iteration, prevents excessive shrinking in first iteration.
#' @param mincumscale double: for types "similarity" and "affine", minimum cumulative scale factor, prevents shrinking below threshold (e.g. 0.7), values above 1 will enforce size increase, should only be used if both landmark configurations are in the same scale.
#' @return landmark array: moving points after the alignment
#' @examples
#' print("pseudotest")
#' @importFrom Rvcg vcgClostKD
#' @export
pts2meshICP <- function(x,
                   # y,
                   maxiter=100,
                   mindist=1e15,
                   minchange=0.0001,
                   type=c("rigid","similarity","affine"),
                   # threads=1,
                   minscale=0.9,
                   mincumscale=NULL) {
  type <- match.arg(type,c("rigid","similarity","affine"))
  cumscalef <- 1
  i=1
  while(i <= maxiter) {
    yclost <- vert2points( vcgClostKD(x, y) )
    distances <- sqrt(rowSums((yclost - x)^2))
    good <- which(distances < mindist)
    trafo <- computeTransform(yclost[good,],x[good,],type=type)
    if(type!="rigid"){
      # calculate scale factor from transformation matrix
      scalef <- sqrt(sum(trafo[1:3,1]^2))
      if (scalef<minscale){
        trafo <- computeTransform(yclost[good,],x[good,], type = "rigid")
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
    xold <- x
    x <- applyTransform(x,trafo)
    eucldist <- sqrt(rowSums((xold - x)^2))
    if(max(eucldist) < minchange){
      i=maxiter
    }
    i = i + 1
  }
  return(x)
}