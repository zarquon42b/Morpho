#' calculate angle between two vectors
#' 
#' calculates unsigned angle between two vectors
#' 
#' 
#' @param x numeric vector (or matrix to be interpreted as vector)
#' @param y numeric vector (or matrix to be interpreted as vector) of same
#' length as \code{x}
#' @return angle between x and y in radians.
#' 
#' @examples
#' 
#' #calculate angle between two centered and
#' # superimposed landmark configuration
#' data(boneData)
#' opa <- rotonto(boneLM[,,1],boneLM[,,2])
#' angle.calc(opa$X, opa$Y)
#' 
#' @export
 angle.calc <- function(x,y)
 {  x <- as.vector(x)/sqrt(sum(x^2))
    y <- as.vector(y)/sqrt(sum(y^2))
    rho <- acos((sum((x-y)^2)-2)/-2)
    return(rho)
 }
