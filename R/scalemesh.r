scalemesh<-function(mesh,size)
{	
	mesh$vb[1:3,]<-mesh$vb[1:3,]*size
	return(mesh)
}
	
