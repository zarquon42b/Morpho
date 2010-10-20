gdif<-function (a3) 
{
    cc <- c.size(addo(a3)/dim(a3)[3])
    x <- sweep(a3, c(1, 2), apply(a3, c(1, 2), mean))
    z <- sqrt(sum((as.vector(x)/cc)^2/dim(a3)[3])^2)
	
    return(z)
}

