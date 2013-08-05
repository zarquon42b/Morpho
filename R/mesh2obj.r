mesh2obj <- function(x,filename="default")
{	x.vb <- cbind("v",t(x$vb[1:3,]))
	x.it <- cbind("f",t(x$it))
	obj <- rbind(x.vb,x.it)
	write.obj(format(obj,scientific=FALSE,trim=TRUE),filename=filename)
}

