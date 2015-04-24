#' calculate Centroid Size for a landmark configuration
#' 
#' calculate Centroid Size for a landmark configuration
#' 
#' 
#' @param x matrix where each row contains coordinates for landmarks
#' @return returns Centroid size
#' 
#' @examples
#' 
#' data(boneData)
#' cSize(boneLM[,,1])
#' 
#' @export
cSize <- function(x)
{	X <- scale(x, scale = FALSE)
	y <- sqrt(sum(as.vector(X)^2))
	return(y)
}
