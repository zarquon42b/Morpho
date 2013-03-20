conv2backf<-function(mesh)
{ 	
	mesh$it<-mesh$it[c(3,2,1),]
        mesh <- adnormals(mesh)
  	return(mesh)
}
