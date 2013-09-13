meshDist <- function(x,...) UseMethod("meshDist")

meshDist.mesh3d <- function(x, mesh2=NULL, distvec=NULL, from=NULL, to=NULL, steps=20, ceiling=FALSE, file="default", imagedim="100x800", uprange=1, ray=FALSE, raytol=50, save=FALSE, plot=TRUE, sign=TRUE, tol=NULL, displace=FALSE, shade=TRUE, method=c("morpho", "vcglib"), ...)
  {
    method=substring(method[1],1L,1L)
    neg=FALSE
    ramp <- blue2green2red(steps-1)
    if (is.null(distvec)) {
        if(!ray) {
            if (method == "v") {
                promesh <- projRead(t(x$vb[1:3,]),mesh2,readnormals=T,sign=T)
            } else {
                promesh <- closemeshKD(x,mesh2,sign=T)
            }
            clost <- promesh$vb[1:3,]
            dists <- promesh$quality
            distsOrig <- dists
            if (!sign)
                dists <- abs(dists)
        } else {
            promesh <- ray2mesh(x,mesh2,tol=raytol,mindist=TRUE)
            clost <- promesh$vb[1:3,]
            dists <- promesh$quality
            distsOrig <- dists
            if (!sign)
              dists <- abs(dists)
        }
    } else {
        clost <- NULL
        dists <- distvec
        distsOrig <- dists
        if (!sign)
            dists <- abs(dists)
    }  
    
    if (is.null(from)) {
        mindist <- min(dists)
        if (sign && mindist < 0 ) {
            from <- quantile(dists,probs=(1-uprange)) 
            neg <- TRUE            
        } else {
            from <- 0
        }
    }
    if (from < 0)
        neg <- TRUE
    if (is.null(to))
        to <- quantile(dists,probs=uprange)    
    if(ceiling)
        to <- ceiling(to)
    
    to <- to+1e-10
    colseq <- seq(from=from,to=to,length.out=steps)
    coldif <- colseq[2]-colseq[1]
    if (neg && sign) {
        negseq <- length(which(colseq<0))
        poseq <- steps-negseq
        maxseq <- max(c(negseq,poseq))
        ramp <- blue2green2red(maxseq*2)
        ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
        distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
        distqual[which(distqual < 1)] <- steps+10
    } else {
        distqual <- ceiling((dists/coldif)+1e-14)
    }
    colorall <- ramp[distqual]
    
    if (!is.null(tol)) {
        if (sign) {
            tol <- c(-tol,tol)
        } else {
            tol <- c(0,tol)
        }
        good <- which(abs(dists) < tol[2])
        colorall[good] <- "#00FF00"
    }   
    
    colfun <- function(x){x <- colorall[x];return(x)}
    x$material$color <- matrix(colfun(x$it),dim(x$it))
    colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
    params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol)
    out <- list(colMesh=x,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
    class(out) <- "meshDist"

    if (plot)
        render(out,output=FALSE,displace=displace,shade=shade,...)
    if (save)
        export(out,file=file,imagedim=imagedim)
    invisible(out)
}
render <- function(x,...) UseMethod("render")
render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,displace=FALSE,shade=TRUE,sign=NULL,...)
  {
    clost <- x$clost
    dists <- x$dists
    distsOrig <- dists
    colorall <- x$cols
    colramp <- x$colramp
    params <- x$params
    distqual <- x$distqual
    if (rgl.cur() !=0)
        rgl.clear()
    
    if (!is.null(from) || !is.null(to) || !is.null(to) || !is.null(uprange) ||  !is.null(tol)  ||  !is.null(sign)) {
        neg=FALSE
        colMesh <- x$colMesh
        if(is.null(steps))
            steps <- x$params$steps
        if(is.null(sign))
          sign <- x$params$sign
        if (!sign) {
            distsOrig <- dists
            dists <- abs(dists)
        }
        if(is.null(ceiling))
            ceiling <- x$params$ceiling
        if(is.null(uprange))
          uprange <- x$params$uprange
                 
        if (is.null(from)) {
            mindist <- min(dists)
            if (sign && mindist < 0 ) {
                from <- quantile(dists,probs=(1-uprange)) 
                neg <- TRUE            
            } else {
                from <- 0
            }             
        }
        if (from < 0)
            neg <- TRUE
        if (is.null(to))
            to <- quantile(dists,probs=uprange)    
        if(ceiling)
            to <- ceiling(to)
        
        to <- to+1e-10
        ramp <- blue2green2red(steps-1)
        colseq <- seq(from=from,to=to,length.out=steps)
        coldif <- colseq[2]-colseq[1]
        if (neg && sign) {
            negseq <- length(which(colseq<0))
            poseq <- steps-negseq
            maxseq <- max(c(negseq,poseq))
            ramp <- blue2green2red(maxseq*2)
            ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
            distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
            distqual[which(distqual < 1)] <- steps+10
        } else {
            distqual <- ceiling((dists/coldif)+1e-14)
        }
        colorall <- ramp[distqual]
        if (!is.null(tol)) {
            if (sign) {
                tol <- c(-tol,tol)
            } else {
                tol <- c(0,tol)
            }
            good <- which(abs(dists) < tol[2])
            colorall[good] <- "#00FF00"
        }
        colfun <- function(x){x <- colorall[x];return(x)}
        colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
        colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
    } else {
        if (is.null(tol))
            tol <- x$params$tol
        
        colramp <- x$colramp
        colMesh <- x$colMesh
    }
    if (shade)
        shade3d(colMesh,specular="black",...)
    if (displace) {
        dismesh <- colMesh
        vl <- dim(colMesh$vb)[2]
        dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
        dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
        dismesh$material$color <- rbind(colorall,colorall,colorall)
        wire3d(dismesh,lit=FALSE)
    }
    diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
    image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
    if (!is.null(tol)) {
        if (sum(abs(tol)) != 0)
            image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col="green",useRaster=TRUE,add=TRUE)
    }
    params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol)
    out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
                                        #out <- list(colMesh=colMesh,colramp=colramp)
    class(out) <- "meshDist"
    invisible(out)
}

export <- function(x,...)UseMethod("export")
export.meshDist <- function(x,file="default",imagedim="100x800",...)
{
    tol <- x$params$tol
    colramp <- x$colramp
    widxheight <- as.integer(strsplit(imagedim,split="x")[[1]])
    mesh2ply(x$colMesh,col=x$cols,filename=file)
    png(filename=paste(file,".png",sep=""),width=widxheight[1],height=widxheight[2])
    diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
    image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
    if (!is.null(tol)) {
        if (sum(abs(tol)) != 0) {
            image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col="green",useRaster=TRUE,add=TRUE)
        }
    }
    dev.off()
}
