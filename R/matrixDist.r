meshDist.matrix <- function(x,mesh2=NULL,distvec=NULL,from=NULL,to=NULL,steps=20,ceiling=FALSE,uprange=1,plot=TRUE,sign=FALSE,tol=NULL,type=c("s","p"),radius=NULL,...)
  {
    x <- list(vb=t(x),it=(1:dim(x)[1]))
    class(x) <- "mesh3d"
    out <- meshDist(x,mesh2=mesh2,distvec=distvec,from=from,to=to,steps=20,ceiling=ceiling,file=file,uprange=uprange ,save=FALSE,plot=FALSE,sign=sign,tol=tol,...)
    class(out) <- "matrixDist"
    render(out,radius=radius,type=type)
    invisible(out)
  }
render.matrixDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,output=FALSE,uprange=NULL,tol=NULL,type=c("s","p"),radius=NULL,...)
  {
    type=type[1]
    dists <- x$dists
    colorall <- x$colorall
    colramp <- x$colramp
    params <- x$params
    distqual <- x$distqual    
    
    if (!is.null(from) || !is.null(to) || !is.null(to) || !is.null(uprange) ||  !is.null(tol))
      {
        
        neg=FALSE
        dists <- x$dists
        sign <- x$params$sign
        colMesh <- x$colMesh
        if(is.null(steps))
          {steps <- x$params$steps
         }
        if(is.null(ceiling))
          {ceiling <- x$params$ceiling
         }
        if(is.null(uprange))
          {uprange <- x$params$uprange
         }
        
        if (is.null(from))
          {
            mindist <- min(dists)
            if (sign && mindist < 0 )
              {
                from <- quantile(dists,probs=(1-uprange)) 
                                        #from <- mindist
                neg <- TRUE            
              }
            else
              {
                from <- 0
              }             
          }
        if (from < 0)
          {
            neg <- TRUE
          }         
        if (is.null(to))
          {to <- quantile(dists,probs=uprange)    
         }
        if(ceiling)
          {to <- ceiling(to)
         }
        to <- to+1e-10
        ramp <- blue2green2red(steps-1)
        colseq <- seq(from=from,to=to,length.out=steps)
        coldif <- colseq[2]-colseq[1]
        if (neg && sign)
          {
            negseq <- length(which(colseq<0))
            poseq <- steps-negseq
            maxseq <- max(c(negseq,poseq))
            ramp <- blue2green2red(maxseq*2)
            ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
            distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
            distqual[which(distqual < 1)] <- steps+10
          }
        else
          {
            distqual <- ceiling((dists/coldif)+1e-14)
          }
        colorall <- ramp[distqual]
         if (!is.null(tol))
           {
             if (sign)
               {
                 tol <- c(-tol,tol)
               }
             else
               {
                 tol <- c(0,tol)
               }
             good <- which(abs(dists) < tol[2])
             colorall[good] <- "#00FF00"
           }
        colfun <- function(x){x <- colorall[x];return(x)}
        colMesh$material$color <- colfun(colMesh$it)
        colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
      }
    else
      {
        if (is.null(tol))
          {
            tol <- x$params$tol
          }
        colramp <- x$colramp
        colMesh <- x$colMesh
      }
    
        if (type == "p")
          {
            points3d(t(colMesh$vb),col=colMesh$material$color,...)
          }
        else
          {
            if (is.null(radius))
              {
                k<-dim(x$colMesh$vb)[2]
                radius <- (cSize(t(x$colMesh$vb))/sqrt(k))*(1/80)
              }
            spheres3d(t(colMesh$vb),col=colMesh$material$color,radius=radius,...)
          }
    diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
    image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
    if (!is.null(tol))
      {
        if (sum(abs(tol)) != 0)
          {
            image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col="green",useRaster=TRUE,add=TRUE)
          }
      }
    params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol)
    out <- list(colMesh=colMesh,dists=dists,cols=colorall,colramp=colramp,params=params,distqual=distqual)
    #out <- list(colMesh=colMesh,colramp=colramp)
     class(out) <- "matrixDist"
    if(output)
      {invisible(out)
     }
  }
