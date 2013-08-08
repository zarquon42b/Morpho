mesh2obj <- function(x,filename=dataname)
{	x.vb <- cbind("v",t(x$vb[1:3,]))
	x.it <- cbind("f",t(x$it))
	obj <- rbind(x.vb,x.it)
        dataname <- deparse(substitute(x))
        filename <- paste(filename,".obj",sep="")
	write.obj(format(obj,scientific=FALSE,trim=TRUE),filename=filename)
}

