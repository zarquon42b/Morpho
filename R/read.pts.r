#' reads pts files
#' 
#' reads Landmark data exported from the software Landmark from
#' http://graphics.idav.ucdavis.edu/research/EvoMorph
#' 
#' 
#' @param file pts file
#' @param na specifies a value that indicates missing values
#' @return
#' \item{matrix }{matrix containing landmark information rownames will be
#' the names given to the landmarks in Landmark}
#' @seealso \code{\link{read.pts}}
#' 
#' @examples
#' 
#' data(nose)
#' write.pts(shortnose.lm, filename="shortnose")
#' data <- read.pts("shortnose.pts")
#' 
#' @export
read.pts <- function(file="x", na=9999)
{	
	pts <- read.table(file,skip=2)
	rownames(pts) <- pts[,1]
	pts <- as.matrix(pts[,2:4])
        nas <- which(pts == na)
        if (length(nas) > 0)
            pts[nas] <- NA
	return(pts)
}
