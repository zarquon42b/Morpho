meshDist <- function(mesh1,mesh2,from=NULL,to=NULL,steps=20,ceiling=FALSE,file="default",imagedim="100x800",uprange=1,ray=FALSE,raytol=50)
  {
    if(!ray)
      {
        dists <- projRead(t(mesh1$vb[1:3,]),mesh2,readnormals=T)$quality
      }
    else
      {
        dists <- ray2mesh(mesh1,mesh2,tol=raytol)$quality
      }
    ramp <- blue2green2red(steps)
    if (is.null(from))
      {from <- 0
     }
    if (is.null(to))
      {to <- quantile(dists,probs=uprange)    
     }
    if(ceiling)
      {to <- ceiling(to)
     }
    colseq <- seq(from=from,to=to,length.out=steps)
    coldif <- colseq[2]-colseq[1]
    distqual <- ceiling((dists/coldif)+1e-14)
    mesh2ply(mesh1,col=ramp[distqual],filename=file)
    widxheight <- as.integer(strsplit(imagedim,split="x")[[1]])
    png(filename=paste(file,".png",sep=""),width=widxheight[1],height=widxheight[2])
    image(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
    dev.off()
    return(list(dists=dists,cols=ramp[distqual]))
  }
