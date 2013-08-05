read.pts <- function(file="x")
{	
	pts <- read.table(file,skip=2)
	rownames(pts) <- pts[,1]
	pts <- as.matrix(pts[,2:4])
	return(pts)
}
