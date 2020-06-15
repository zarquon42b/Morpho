#' check for NA values in a matrix (of landmarks)
#'
#' check for NA values in a matrix (of landmarks)
#'
#' @param x matrix containing landmarks
#' @return returns a vector with missin landmarks and a vector of length=0 if none are missing
#' @export
checkNA <- function(x) {
    chk <- rowSums(x)
    ignore <- which(is.na(chk))
    return(ignore)
}
        
