estShape<-function(shapelist,procrot,n,ref)
{	
	wcalc<-proc.weight(procrot,n,ref)
	dimvb<-dim(shapelist[[1]]$vb)
	
	vb<-matrix(0,3,dimvb[2])
	
	for (i in 1:n)
		{vb<-vb+shapelist[[wcalc$data$nr[i]]]$vb[1:3,]*wcalc$data$weight[i]
		}
	
	lms<-procrot[,,wcalc$data$nr]
	lm.est<-matrix(0,dim(procrot)[1],3)
	for (i in 1:n)
		{lm.est<-lm.est+lms[,,i]*wcalc$data$weight[i]
		}
	mesh<-shapelist[[1]]
	mesh$vb[1:3,]<-vb	
	return(list(mesh=mesh,lm.est=lm.est))
	#weights<-wcalc$data[,4]
	
}
