deformGrid<-function(matrix,tarmatrix,ngrid=10,lwd=1,showaxis=c(1,2,3),both=T,lines=TRUE,lcol=1,add=FALSE,col1=2,col2=3,type=c("s","p"))
{

  type <- type[1]
  out3d <- spheres3d
  if (type == "p")
    {out3d <- points3d
   }
  if (!add)
    {
      open3d()
    }
  k<-dim(matrix)[1]
  sz <- (cSize(matrix)/sqrt(k))*(1/80)
  out3d(matrix,col=col1,radius=sz)
  if(both)
    {out3d(tarmatrix,col=col2,radius=sz)
     if (lines)
       {
         linemesh <- list()
         linemesh$vb <- t(cbind(rbind(matrix,tarmatrix),1))
         linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
         class(linemesh) <- "mesh3d"
         wire3d(linemesh,lwd=1.5,col=lcol,lit=FALSE)
                                        #               for (i in 1:k)
                                        #lines3d(rbind(matrix[i,],tarmatrix[i,]),lwd=1.5)
       }
   }
  x2<-x1<-x3<-c(0:(ngrid-1))/ngrid;x0<-as.matrix(expand.grid(x1,x2,x3))
  cent.mat<-apply(matrix,2,scale,scale=F)
  mean.mat<-apply(matrix,2,mean)
  
  if (ngrid > 1)
    {
      xrange<-diff(range(matrix[,1]))
      yrange<-diff(range(matrix[,2]))
      zrange<-diff(range(matrix[,3]))
      xrange1<-diff(range(tarmatrix[,1]))
      yrange1<-diff(range(tarmatrix[,2]))
      zrange1<-diff(range(tarmatrix[,3]))
      maxi<-max(c(xrange,yrange,zrange,xrange1,yrange1,zrange1))
      maxi<-maxi+0.02*maxi
      x0<-maxi*x0
      x0<-apply(x0,2,scale,scale=FALSE)
      space<-eigen(crossprod(cent.mat))$vectors
      x0<-t(t(x0%*%space)+mean.mat)
      x0<-tps3d(x0,matrix,tarmatrix)
      
      
      it <- NULL
      if(1 %in% showaxis)	
	{for (j in 0:(ngrid-1))
           {for (i in 0:(ngrid-1))
              {a<-(c((j*ngrid+1):(j*ngrid+ngrid)))*ngrid+i-(ngrid-1)
               it <- cbind(it,a)
               lines3d(x0[a,],smooth=T,col=1,lwd=lwd)
             }
          }
       }
      if(2 %in% showaxis)	
        {for (i in c(0:(ngrid^2-1))*ngrid)
           {lines3d(x0[(i+1):(i+ngrid),],smooth=T,col=1,lwd=lwd)}
       }	
      if(3 %in% showaxis)	
        {for (i in c(0:(ngrid^2-1)))
           {a<-c(1:ngrid)*ngrid^2-i
            lines3d(x0[a,],smooth=T,col=1,lwd=lwd)}
       }
    }
}

