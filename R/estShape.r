estShape<-function(shapelist,procrot,n,ref,softproc=NULL,size=NULL)
{	
	wcalc<-proc.weight(procrot,n,ref)
	dimvb<-dim(shapelist[[1]]$vb)
	lm.est<-NULL
	vb<-matrix(0,3,dimvb[2])
	meandev<-NULL
	for (i in 1:n)
		{vb<-vb+shapelist[[wcalc$data$nr[i]]]$vb[1:3,]*wcalc$data$weight[i]
		print(names(shapelist)[wcalc$data$nr[i]])
		}
	
	
	if (! is.null(softproc))
	{
	if (is.null(size))
	stop("please enter size vector")
	
	lms<-softproc[,,wcalc$data$nr]
	lm.est<-matrix(0,dim(softproc)[1],3)
	if (i > 1)
	{	
	for (i in 1:n)
		{lm.est<-lm.est+lms[,,i]*wcalc$data$weight[i]
		}
	}
	else 
	{lm.est<-lms
	}
	meandev<-mean(sqrt(diag(tcrossprod((lm.est-softproc[,,ref])*size[ref]))))
	
	}
	mesh<-shapelist[[1]]
	mesh$vb[1:3,]<-vb
	mesh<-adnormals(mesh)	
	return(list(mesh=mesh,lm.est=lm.est,weights=wcalc$data,meandev=meandev))
	
	
}
