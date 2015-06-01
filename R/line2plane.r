#' get intersection between a line and a plane
#'
#' get intersection between a line and a plane
#' @param ptLine vector of length 3: point on line
#' @param ptDir  vector of length 3: direction vector of line
#' @param planeNorm vector of length 3: plane normal vector
#' @param planePt vector of length 3: point on plane
#' @return hit point
#' @note in case you only have three points on a plane (named \code{pt1, pt2, pt3} you can get the plane's normal by calling \code{crossProduct(pt1-pt2,pt1-pt3)}.
#' @export
line2plane <- function(ptLine,ptDir, planePt, planeNorm) {
    d <- crossprod(planeNorm,planePt)
    t <- (d-crossprod(planeNorm,ptLine))/crossprod(planeNorm,ptDir)
    out <- ptLine+t*ptDir
    if (!length(out))
        stop("plane and line are collinear")
    return(out)
}
