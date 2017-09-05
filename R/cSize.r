#' calculate Centroid Size for a landmark configuration
#' 
#' calculate Centroid Size for a landmark configuration
#' 
#' 
#' @param x k x 3 matrix containing landmark coordinates or mesh of class "mesh3d" 
#' @return returns Centroid size
#' 
#' @examples
#' 
#' data(boneData)
#' cSize(boneLM[,,1])
#' 
#' @export
cSize <- function(x){
    if(inherits(x,"mesh3d"))
        x <- vert2points(x)
    X <- scale(x, scale = FALSE)
    y <- sqrt(sum(as.vector(X)^2))
    return(y)
}
