#' create a list with empty entries to be used as missingList in slider3d
#'
#' create a list with empty entries to be used as missingList in slider3d
#'
#' @param x length of the list to be created
#' @return returns a list of length \code{x} filled with numerics of length zero.
#' @examples
#' ## Assume in a sample of 10, the 9th individual has (semi-)landmarks 10:50
#' #   hanging in thin air (e.g. estimated using fixLMtps)
#' #   while the others are complete.
#' ## create empty list
#' missingList <- createMissingList(10)
#' missingList[[9]] <- 10:50
#' @seealso \code{\link{fixLMtps},\link{fixLMmirror}, \link{slider3d}} 
#' @export
createMissingList <- function(x) {
    ml <- lapply(1:x,function(x) x <- numeric(0))
    return(ml)
}
