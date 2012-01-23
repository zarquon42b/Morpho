mc.rot<-function(dat.array,tol=1e-5)
{
	
	x<-dat.array
	p1<-1e10	
	n<-dim(dat.array)[3]
	k<-dim(dat.array)[1]
	m<-dim(dat.array)[2]
	x1<-gdif(dat.array)
	arr.list<-list(0)	
	for ( i in 1:n)
		{arr.list[[i]]<-list(x[,,i],1)
		}
	mshape<-x[,,1]	
	arr.list<-mclapply(arr.list,function(x){x[[1]]<-apply(x[[1]],2,scale,scale=F);return(list(x[[1]],x[[2]]))})
	while (p1 > tol)
		{
		mshape_old<-mshape

	### rotation of all configs on current consensus ###		
		arr.list<-mclapply(arr.list,function(x){x[[1]]<-rot.proc(x[[1]],x=mshape,scale=F);return(list(x[[1]],x[[2]]))})
		
					
			for( i in 1:n)
			{x[,,i]<-arr.list[[i]][[1]]
			}
		
		x2<-gdif(x)
		p1<-abs(x1-x2)
		x1<-x2
	
		}
	return(x)
}
