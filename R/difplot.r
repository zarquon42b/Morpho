difplot<-function(refshape,targetshape,ID="default",away.fac=1.1,color=4,lwd=1,lcol=2,sz=0,rgl.new=TRUE)
{   
    if ( rgl.new == TRUE)
       {
         rgl.open()
         rgl.bg(color = "white")
       }
    A<-refshape
    k<-dim(A)[1]
    m<-dim(A)[2]
    texx<-apply(targetshape,2,max)
    #texm<-apply(refshape,2,mean)
    texx<-texx*away.fac
    sds<-0
     lim<-max(abs(refshape))
    
    
      if(sz==0)
      {sz <- centroid.size(refshape)/sqrt(k)/30}

      #if (!spheres)
      spheres3d(refshape,  col = color,radius=sz)
      #else {spheres3d(refshape, radius = sz, col = color)}
      #plot3d(refshape,box=F,axes=F,type=type,xlim=c(-lim,lim),ylim=c(-lim,lim),zlim=c(-lim,lim),col=color,xlab="",ylab="",zlab="",radius=sz)
    for (j in 1:k)
          {lines3d(rbind(refshape[j,],targetshape[j,]),col=lcol,lwd=lwd)}
      text3d(texx,text=ID,cex=100,font=5,col=lcol) 
    
}
