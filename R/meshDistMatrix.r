#' @rdname meshDist
#' @method meshDist matrix
#' @export
meshDist.matrix <- function(x,mesh2=NULL,distvec=NULL,from=NULL,to=NULL,steps=20,ceiling=FALSE, rampcolors=colorRamps::blue2green2red(steps-1),NAcol="white", uprange=1,plot=TRUE,sign=TRUE,tol=NULL,type=c("s","p"),radius=NULL,displace=FALSE,add=FALSE,scaleramp=FALSE,...)
    {
        x <- list(vb=t(x),it=matrix(1:dim(x)[1]),1,dim(x)[1])
        class(x) <- "mesh3d"
        out <- meshDist(x,mesh2=mesh2,distvec=distvec,from=from,to=to,steps=20,ceiling=ceiling,file=file,uprange=uprange ,save=FALSE,plot=FALSE,sign=sign,tol=tol,rampcolors = rampcolors,displace=FALSE,NAcol = NAcol,scaleramp=scaleramp,...)
        class(out) <- "matrixDist"
        render(out,radius=radius,type=type,displace=displace,add=add)
        invisible(out)
    }
#' @rdname render
#' @method render matrixDist
#' @export
render.matrixDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,type=c("s","p"),radius=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,sign=NULL,add=FALSE,scaleramp=FALSE,...) {
    if (!add) {
        if (rgl.cur() !=0)
            rgl.clear()
        clost=x$clost
    }
    type=type[1]
    dists <- x$dists
    distsOrig <- dists
    colorall <- x$cols
    colramp <- x$colramp
    params <- x$params
    distqual <- x$distqual    
    
    if (!is.null(from) || !is.null(to) || !is.null(to) || !is.null(uprange) ||  !is.null(tol) || !is.null(sign) || !is.null(rampcolors) || !is.null(NAcol) || !is.null(scaleramp)) {
        neg=FALSE
        dists <- x$dists
        distsOrig <- dists
        colMesh <- x$colMesh
        if(is.null(steps))
            steps <- x$params$steps
        if(is.null(sign))
            sign <- x$params$sign
        if (is.null(rampcolors))
            rampcolors <- x$params$rampcolors
        if (is.null(NAcol))
            NAcol <- x$params$NAcol
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
        if (is.null(scaleramp))
            scaleramp <- x$scaleramp
        if (from < 0)
            neg <- TRUE
        if (is.null(to))
            to <- quantile(dists,probs=uprange)    
        if(ceiling)
            to <- ceiling(to)
        to <- to+1e-10
        ramp <- colorRampPalette(rampcolors)(steps-1)
        colseq <- seq(from=from,to=to,length.out=steps)
        coldif <- colseq[2]-colseq[1]
        if (neg && sign) {
            negseq <- length(which(colseq<0))
            poseq <- steps-negseq
            maxseq <- max(c(negseq,poseq))
            if (scaleramp) {
                ramp <- colorRampPalette(rampcolors)(maxseq*2)
                ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
                
            }
            else
                ramp <- colorRampPalette(rampcolors)(steps-1)
            distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
            distqual[which(distqual < 1)] <- steps+10
        } else if (from > 0) {
              distqual <- ceiling(((dists-from)/coldif)+1e-14)
          } else {
                distqual <- ceiling((dists/coldif)+1e-14)
            }
        distqual[which(distqual < 1)] <- steps+10
        colorall <- ramp[distqual]
        if (!is.null(tol)) {
            if (sign)
                tol <- c(-tol,tol)
            else
                tol <- c(0,tol)
            good <- which(abs(dists) < tol[2])
            colorall[good] <- "#00FF00"
        }
        colfun <- function(x){x <- colorall[x];return(x)}
        colMesh$material$color <- colfun(colMesh$it)
        colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
        
        colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
    } else {
        if (is.null(tol))
            tol <- x$params$tol
        colramp <- x$colramp
        colMesh <- x$colMesh
    }
    
    if (type == "p") {
        points3d(t(colMesh$vb),col=colMesh$material$color,...)
    } else {
        if (is.null(radius)) {
            k <- dim(x$colMesh$vb)[2]
            radius <- (cSize(t(x$colMesh$vb))/sqrt(k))*(1/80)
        }
        spheres3d(t(colMesh$vb),col=colMesh$material$color,radius=radius,...)
    }
    if (displace) {
        dismesh <- colMesh
        vl <- dim(colMesh$vb)[2]
        dismesh$vb <- cbind(rbind(colMesh$vb,1),rbind(clost,1))
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
    out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,NAcol=NAcol)
                                        #out <- list(colMesh=colMesh,colramp=colramp)
    class(out) <- "matrixDist"
    invisible(out)
}
