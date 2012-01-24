difplot.lm2D<-function(refshape,targetshape,color=4,lwd=1,lcol=2,main=main)
{   
    A<-refshape
    k<-dim(A)[1]
    m<-dim(A)[2]
    
    sds<-0
     #lim<-max(abs(refshape))
    
    
      
      #sz <- (cSize(refshape)/sqrt(k))*(1/80)
	
      #if (!spheres)
      plot(refshape,  col = color,main=main,asp=1,axes=FALSE,xlab="",ylab="")
      #else {spheres3d(refshape, radius = sz, col = color)}
      #plot3d(refshape,box=F,axes=F,type=type,xlim=c(-lim,lim),ylim=c(-lim,lim),zlim=c(-lim,lim),col=color,xlab="",ylab="",zlab="",radius=sz)
    for (j in 1:k)
          {lines(rbind(refshape[j,],targetshape[j,]),col=lcol,lwd=lwd)}
      text(refshape,labels=paste("",c(1:k),sep=""),col=lcol,pos=2) 
    
}
