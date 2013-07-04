pcaplot3d <- function (x,...) UseMethod("pcaplot3d")
pcaplot3d.symproc<-function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,sym=TRUE,...) 
{   
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
  if (is.vector(PCs))
      PCs <- matrix(PCs,length(PCs),1)
  
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

  for (i in 1:length(pcshow))
    {pc<-refshape+matrix(PCs[,pcshow[i]]*mag[i]*sds[pcshow[i]],k,3)

     linemesh <- list()
     linemesh$vb <- t(cbind(rbind(refshape,pc),1))
     linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
     class(linemesh) <- "mesh3d"
     wire3d(linemesh,lwd=lwd,lit=F,col=rainb[i])
   }
  spheres3d(refshape,  col = color,radius=sz)
}
pcaplot3d.nosymproc<-function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,...)
{   
                                        #rgl.open()
                                        #rgl.bg(color = "white")
  refshape <- x$mshape
  PCs <- x$PCs
  Scores <- x$PCscores
  A<-refshape
  k<-dim(A)[1]
  m<-dim(A)[2]
  if (is.vector(PCs))
      PCs <- matrix(PCs,length(PCs),1)
 
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

  for (i in 1:length(pcshow))
    {pc<-refshape+matrix(PCs[,pcshow[i]]*mag[i]*sds[pcshow[i]],k,3)

     linemesh <- list()
     linemesh$vb <- t(cbind(rbind(refshape,pc),1))
     linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
     class(linemesh) <- "mesh3d"
     wire3d(linemesh,lwd=lwd,lit=F,col=rainb[i])
   }
  spheres3d(refshape,  col = color,radius=sz)
}
