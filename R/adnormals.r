adnormals<- function(x) 
{	#### this is basically the function addNormals from the rgl-package with an option for checking degenerated faces ####
	

	veclen <- function(v) sqrt(sum(v^2))
	normalize <- function(v) v/veclen(v)
  v <- x$vb
  
  # Make sure v is homogeneous with unit w
  if (nrow(v) == 3) v <- rbind(v, 1)
  else v <- t( t(v)/v[4,] )
  
  normals <- v*0
  v <- v[1:3,]
  
  
    it <- x$it
	storage.mode(v)<-"double"
	storage.mode(normals)<-"double"
	storage.mode(it)<-"integer"
    out<-.Fortran("adnormals",v,ncol(v),it,ncol(it),nrow(it),normals=normals)$normals
	

  normals <- out#/out[4,]

  x$normals <- normals
  x
}


