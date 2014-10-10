#' create a perfectly symmetric version of landmarks
#'
#' create a perfectly symmetric version of landmarks
#'
#' @param x k x m matrix with rows containing landmark coordinates
#' @param pairedLM A X x 2 matrix containing the indices (rownumbers) of the
#' paired LM. E.g. the left column contains the lefthand landmarks, while the
#' right side contains the corresponding right hand landmarks.
#' @return a symmetrized version of \code{x}
#' @details the landmarks are reflected and relabled according to
#' \code{pairedLM} and then rotated and translated onto \code{x}.
#' Both configurations are then averaged to obtain a perfectly symmetric one.
#' @references
#' Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape analysis of symmetric
#' structures: quantifying variation among individuals and asymmetry. Evolution
#' 56(10):1909-1920.
#' @examples
#' data(boneData)
#' left <- c(4,6,8)
#' right <- c(3,5,7)
#' pairedLM <- cbind(left,right)
#' symx <- symmetrize(boneLM[,,2],pairedLM)
#' \dontrun{
#' deformGrid3d(symx,boneLM[,,2])
#' }
#' @export
symmetrize <- function(x, pairedLM) {
    xmir <- mirror(x,icpiter=0)
    xmir[c(pairedLM),] <- xmir[c(pairedLM[,2:1]),]
    xrot <- rotonto(x,xmir)$yrot
    xsym <- (x+xrot)/2
    return(xsym)
}
