deform.grid<-function(matrix,tarmatrix,ngrid,lwd=1,showaxis=c(1,2,3))
{	spheres3d(matrix,col=2)
	x2<-x1<-x3<-c(0:(ngrid-1))/ngrid;x0<-as.matrix(expand.grid(x1,x2,x3))
	
	xrange<-diff(range(matrix[,1]))
	yrange<-diff(range(matrix[,2]))
	zrange<-diff(range(matrix[,3]))
	maxi<-max(c(xrange,yrange,zrange))
	x0<-maxi*x0
	x0<-apply(x0,2,scale,scale=FALSE)
	space<-eigen(crossprod(matrix))$vectors
	x0<-x0%*%space
	x0<-tps3d(x0,matrix,tarmatrix)
	
	if(1 %in% showaxis)	
	{for (j in 0:(ngrid-1))
		{for (i in 0:(ngrid-1))
			{a<-(c((j*ngrid+1):(j*ngrid+ngrid)))*ngrid+i-(ngrid-1);lines3d(x0[a,],smooth=T,col=1,lwd=lwd)
			}
		}
	}
	if(2 %in% showaxis)	
		{for (i in c(0:(ngrid^2-1))*ngrid){lines3d(x0[(i+1):(i+ngrid),],smooth=T,col=1,lwd=lwd)}
		}	
	if(3 %in% showaxis)	
		{for (i in c(0:(ngrid^2-1))){a<-c(1:ngrid)*ngrid^2-i;lines3d(x0[a,],smooth=T,col=1,lwd=lwd)}
		}
}

