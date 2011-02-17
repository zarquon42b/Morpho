vecx<-function(x)
{  dims<-dim(x)
	n<-dims[3]
	k<-dims[1]
	m<-dims[2]
	vecs<-matrix(0,n,k*m)
	for(i in 1:n)
		{vecs[i,]<-as.vector(x[,,i])
		}
	#vecs<-apply(vecs,2,scale,scale=F)
	return(vecs)
}
	
