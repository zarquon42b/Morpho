#' @rdname tangentPlane
#' @export
crossProduct <- function(x,y)
{	
	out <- c(0,0,0)
 	out[1] <- x[2]*y[3]-x[3]*y[2]
	out[2] <- x[3]*y[1]-x[1]*y[3]
	out[3] <- x[1]*y[2]-x[2]*y[1]
        out <- out/norm(out,"2")
 return(out)
 }

