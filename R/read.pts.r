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
