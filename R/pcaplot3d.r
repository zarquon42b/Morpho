pcaplot3d <- function(x,pcshow=c(1,2,3),mag=3,color=3,lwd=1, ...) UseMethod("pcaplot3d")
pcaplot3d.symproc<-function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,sym=TRUE)
{   
                                        #rgl.open()
                                        #rgl.bg(color = "white")
  refshape <- x$mshape
  if (sym)
    {
      PCs <- x$PCsym
      Scores <- x$PCscore_sym
    }
  else
    {
      PCs <- x$PCasym
      Scores <- x$PCscore_asym
    }
  A<-refshape
  k<-dim(A)[1]
  m<-dim(A)[2]
  npc<-dim(PCs)[2]
  lpc<-length(pcshow)
  rainb<-rainbow(lpc)
  sds<-0
  if (length(mag)==1)
    {mag<-c(rep(mag,lpc))}
    lim<-max(abs(refshape))
    
          for (i in 1:npc)
          {sds[i]<-sd(Scores[,i])}

      sz <- cSize(refshape)/sqrt(k)*(1/80)

      #if (!spheres)
      spheres3d(refshape,  col = color,radius=sz)
      #else {spheres3d(refshape, radius = sz, col = color)}
      #plot3d(refshape,box=F,axes=F,type=type,xlim=c(-lim,lim),ylim=c(-lim,lim),zlim=c(-lim,lim),col=color,xlab="",ylab="",zlab="",radius=sz)
  for (i in 1:length(pcshow))
      {pc<-refshape+matrix(PCs[,pcshow[i]]*mag[i]*sds[pcshow[i]],k,3)

          for (j in 1:k)
          {lines3d(rbind(refshape[j,],pc[j,]),col=rainb[i],lwd=lwd)}
      }
}
pcaplot3d.nosymproc<-function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1)
{   
    #rgl.open()
    #rgl.bg(color = "white")
   refshape <- x$mshape
   PCs <- x$PCs
   Scores <- x$PCscores
    A<-refshape
    k<-dim(A)[1]
    m<-dim(A)[2]
    npc<-dim(PCs)[2]
    lpc<-length(pcshow)
    rainb<-rainbow(lpc)
    sds<-0
    if (length(mag)==1)
      {mag<-c(rep(mag,lpc))}
    lim<-max(abs(refshape))
    
          for (i in 1:npc)
          {sds[i]<-sd(Scores[,i])}

      sz <- cSize(refshape)/sqrt(k)*(1/80)

      #if (!spheres)
      spheres3d(refshape,  col = color,radius=sz)
      #else {spheres3d(refshape, radius = sz, col = color)}
      #plot3d(refshape,box=F,axes=F,type=type,xlim=c(-lim,lim),ylim=c(-lim,lim),zlim=c(-lim,lim),col=color,xlab="",ylab="",zlab="",radius=sz)
  for (i in 1:length(pcshow))
      {pc<-refshape+matrix(PCs[,pcshow[i]]*mag[i]*sds[pcshow[i]],k,3)

          for (j in 1:k)
          {lines3d(rbind(refshape[j,],pc[j,]),col=rainb[i],lwd=lwd)}
      }
}
