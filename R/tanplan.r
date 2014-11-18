#' calculate the orthogonal complement of a 3D-vector
#' 
#' calculate the orthogonal complement of a 3D-vector or the 3D-crossproduct,
#' finding an orthogonal vector to a plane in 3D.
#' 
#' 
#' @title calculate the orthogonal complement of a 3D-vector
#' @param x vector of length 3.
#' @param y vector of length 3.
#' @return tangentPlane:
#' 
#' crossProduct: returns a vector of length 3.
#' \item{y }{vector orthogonal to x}
#' \item{z }{vector orthogonal to x and y}
#' @author Stefan Schlager
#' @examples
#' 
#' require(rgl)
#' 
#' x <- c(1,0,0)
#' y <- c(0,1,0)
#' 
#' #example tangentPlane
#' z <- tangentPlane(x)
#' #visualize result
#' \dontrun{
#' lines3d(rbind(0, x), col=2, lwd=2)
#' ## show complement
#' lines3d(rbind(z$y, 0, z$z), col=3, lwd=2)
#' }
#' # example crossProduct
#' z <- crossProduct(x, y)
#' # show x and y
#' \dontrun{
#' lines3d(rbind(x, 0, y), col=2, lwd=2)
#' # show z
#' lines3d(rbind(0, z), col=3, lwd=2)
#' }
#' @rdname tangentPlane
#' @export
tangentPlane <- function(x)
{		if (sum(x^2)==0)
		{stop(cat("zero vector has no orthogonal subspace"))}
		
		
		
		if (0 %in% x)
			{
			y <- c(0,0,0)	
			y[which(x==0)] <- 1
                        y <- y/sqrt(sum(y^2))
			}
		else 
			{
			y <- c(1,1,-(x[1]+x[2])/x[3])
			y <- y/sqrt(sum(y^2))
			}
		z <- crossProduct(y,x)
		z <- z/sqrt(sum(z^2))
	return(list(z=z,y=y))
}
			
		
