#' convert an 3D array into a matrix
#' 
#' converts a 3D-array (e.g. containing landmark coordinates) into a matrix,
#' one row per specimen.
#' 
#' 
#' @param x array
#' 
#' @param byrow logical: if TRUE, the resulting vector for each specimen will
#' be x1,y1,z1,x2,y2,z2,..., and x1,x2,...,y1,y2,...,z1,z2,... otherwise
#' (default)
#' 
#' @return returns a matrix with one row per specimen
#' @author Stefan Schlager
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' library(shapes)
#' data <- vecx(gorf.dat)
#' 
#' @export vecx
vecx <- function(x,byrow=FALSE)
{  dims <- dim(x)
	n <- dims[3]
	k <- dims[1]
	m <- dims[2]
   names <- dimnames(x)[[3]]
   
	vecs <- matrix(0,n,k*m)
	for(i in 1:n)
          {
            if (byrow)
              {
                vecs[i,] <- as.vector(t(x[,,i]))
              }
            else
              {vecs[i,] <- as.vector(x[,,i])
             }
          }
   if (!is.null(names))
     {
       rownames(vecs) <- names
     }
	#vecs <- apply(vecs,2,scale,scale=F)
	return(vecs)
}
	
