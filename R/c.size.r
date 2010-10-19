c.size<-function(x)
{	X <- apply(x, 2, scale, scale = F)
	y<- sqrt(sum(as.vector(X)^2))
	return(y)
}
