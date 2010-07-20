read.pts<-function(file="x")
{	
	pts<-read.table(file,skip=2)
	rownames(pts)<-pts[,1]
	pts<-pts[,2:4]
	return(pts)
}
