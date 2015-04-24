#' convert an 3D array into a matrix and back
#' 
#' converts a 3D-array (e.g. containing landmark coordinates) into a matrix,
#' one row per specimen or reverse this.
#' 
#' 
#' @param x array or matrix
#' 
#' @param byrow logical: if TRUE, the resulting vector for each specimen will
#' be \code{x1,y1,z1,x2,y2,z2,...,} and \code{x1,x2,...,y1,y2,...,z1,z2,...} otherwise
#' (default). The same is for reverting the process: if the matrix contains the coordinates as  rows like: \code{x1,y1,z1,x2,y2,z2,...} set \code{byrow=TRUE}
#' @param revert revert the process and convert a matrix with vectorized landmarks back into an array.
#' @param lmdim number of columns for reverting
#' 
#' @return returns a matrix with one row per specimen
#' @author Stefan Schlager
#' 
#' @examples
#' 
#' library(shapes)
#' data <- vecx(gorf.dat) 
#' #revert the procedure
#' gdat.restored <- vecx(data,revert=TRUE,lmdim=2)
#' range(gdat.restored-gorf.dat)
#' @export
vecx <- function(x, byrow=FALSE, revert=FALSE, lmdim) {
    dims <- dim(x)
    if (!revert) {
        n <- dims[3]
        k <- dims[1]
        m <- dims[2]
    } else {
        n <- dims[1]
        k <- ncol(x)/lmdim
        m <- lmdim
    }
    if (!revert)
        names <- dimnames(x)[[3]]
    else
        names <- rownames(x)
    if (!revert) {
        vecs <- matrix(0,n,k*m)
                                        #for(i in 1:n) {
        if (byrow) {
            vecs <- matrix(aperm(x,c(3,2,1)),n,k*m)
        } else {
            vecs <- matrix(aperm(x,c(3,1,2)),n,k*m)
        }
    } else {
        x <- t(x)
        vecs <- array(x, dim=c(k,m,n))
        if (byrow) {
            tmp <- array(x,dim=c(m,k,n))
            for (i in 1:n)
                vecs[,,i] <- t(tmp[,,i])
        }
    }
    if (!is.null(names) && !revert) {
        rownames(vecs) <- names
    }
    return(vecs)
}

