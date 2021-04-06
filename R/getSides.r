#' try to identify bilateral landmarks and sort them by side
#'
#' try to identify bilateral landmarks and sort them by side
#' @param x matrix containing landmarks (see details)
#' @param tol maximal distance allowed between original and mirrored set.
#' @param ... additional arguments passed to \code{\link{mirror}}.
#' @return returns a list containing
#' \item{side1}{integer vector containing indices of landmarks on one side}
#' \item{side2}{integer vector containing indices of landmarks on the other side}
#' \item{unilat}{integer vector containing indices unilateral landmarks}
#' @details This function mirrors the landmark set and aligns it to the original. Then it tries to find pairs. If you have a sample, run a Procrustes registration first (without scaling to unit centroid size, or you later have to adapt \code{tol} - see examples) and then use the mean as it is usually more symmetrical.
#' @examples
#' data(boneData)
#' proc <- procSym(boneLM,CSinit=FALSE)
#' mysides <- getSides(proc$mshape)
#' if (interactive()){
#' #visualize bilateral landmarks
#' deformGrid3d(boneLM[mysides$side1,,1],boneLM[mysides$side2,,1])
#' ## visualize unilateral landmarks
#' rgl::spheres3d(boneLM[mysides$unilat,,1],radius=0.5)
#' }
#' @importFrom Rvcg vcgKDtree
#' @export
getSides <- function(x,tol=3,...) {
    xmir <- mirror(x,...)
    clost <- vcgKDtree(x,xmir,k=1)
    myindex <- (1:nrow(x))
    pairs <- which(clost$distance < tol)
    side1 <- myindex[pairs]
    side2 <- clost$index[pairs]
    unilat <- which(side1==side2)
    if (length(unilat)) {
        side1 <- side1[-unilat]
        side2 <- side2[-unilat]
    }
    meanplane <- unique(x[side1,]+x[side2,])/2
    mycut <- cutSpace(pointcloud=x[side1,],v1=meanplane[1,],v2=meanplane[2,],v3=meanplane[3,])
    left <- side1[which(mycut)]
    right <- side2[which(mycut)]
    unilat <- myindex[unilat]

    return(list(side1=left,side2=right,unilat=unilat))
}



