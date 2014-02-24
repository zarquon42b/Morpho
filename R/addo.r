addo <- function(arr) {
    if (!is.numeric(arr) || length(dim(arr)) != 3)
        stop("please provide 3D numeric array")
    out <- .Call("addo",arr)
    return(out)
}
#' calculate mean of an array
#'
#' calculate mean of a 3D-array (e.g. containing landmarks) (fast) using the Armadillo C++ Backend
#'
#' @param arr \code{k x m x n} dimensional numeric array
#' @return matrix of dimensions \code{k x m}.
#' @note this is the same as \code{apply(arr, 1:2, mean)}, only faster for large configurations.
#' @examples
#' data(boneData)
#' proc <- ProcGPA(boneLM, silent = TRUE)
#' mshape <- arrMean3(proc$rotated)
#' @export
arrMean3 <- function(arr) {
    if (!is.numeric(arr) || length(dim(arr)) != 3)
        stop("please provide 3D numeric array")
    out <- .Call("arrMean3",arr)
    return(out)
}
