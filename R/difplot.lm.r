difplot.lm<-function(refshape,targetshape,color=4,lwd=1,lcol=2,rgl.new=TRUE)
{   
    if ( rgl.new == TRUE)
       {
         open3d()
       }
    A<-refshape
    k<-dim(A)[1]
    m<-dim(A)[2]
    #texx<-apply(targetshape,2,max)
    #texm<-apply(refshape,2,mean)
    #texx<-texx*away.fac
    sds<-0
     lim<-max(abs(refshape))
    
    
      
      sz <- (cSize(refshape)/sqrt(k))*(1/80)
	
      #if (!spheres)
      spheres3d(refshape,  col = color,radius=sz)
      #else {spheres3d(refshape, radius = sz, col = color)}
      #plot3d(refshape,box=F,axes=F,type=type,xlim=c(-lim,lim),ylim=c(-lim,lim),zlim=c(-lim,lim),col=color,xlab="",ylab="",zlab="",radius=sz
    linemesh <- list()
    linemesh$vb <- t(cbind(rbind(refshape,targetshape),1))
    linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
    class(linemesh) <- "mesh3d"
    wire3d(linemesh,col=lcol,lwd=lwd,lit=F)
       
    text3d(refshape,texts=paste("",c(1:k),sep=""),cex=1,col=lcol,adj=1.2) 
    
}
