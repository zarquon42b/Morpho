#' visualise differences between two superimposed sets of 3D landmarks
#' 
#' visualise differences between two superimposed sets of 3D landmarks by
#' deforming a cubic grid based on a thin-plate spline interpolation
#' 
#' 
#' @param matrix reference matrix containing 3D landmark coordinates.
#' @param tarmatrix target matrix containing 3D landmark coordinates.
#' @param ngrid number of grid lines to be plotted; ngrid=0 suppresses grid
#' creation.
#' @param lwd width of lines connecting landmarks.
#' @param showaxis integer (vector): which dimensions of the grid to be
#' plotted. Options are combinations of 1,2 and 3.
#' @param both logical: if FALSE, only "matrix" will be plotted.
#' @param lines logical: if TRUE, lines between landmarks will be plotted.
#' @param lcol color of lines
#' @param add logical: add to existing rgl window.
#' @param col1 color of "matrix"
#' @param col2 color of "tarmat"
#' @param type "s" renders landmarks as spheres; "p" as points - much faster
#' for very large pointclouds.
#' @param size control size/radius of points/spheres
#' @author Stefan Schlager
#' @seealso \code{\link{tps3d}}
#' 
#' @examples
#' \dontrun{
#' data(nose)
#' deformGrid3d(shortnose.lm,longnose.lm,ngrid=10)
#' }
#' @export
deformGrid3d <- function(matrix,tarmatrix,ngrid=0,lwd=1,showaxis=c(1, 2), both=T,lines=TRUE,lcol=1,add=FALSE,col1=2,col2=3,type=c("s","p"),size=NULL)
{
    type <- type[1]
    if (dim(matrix)[1] > 1000 && type =="s" && size > 0) {
        answer <- readline("You have a lot of landmarks\n Render them as points (faster)? (yes/NO)\n")
        if (! substr(answer,1L,1L) %in% c("n","N"))
            type <- "p"
    }
    out3d <- spheres3d
    if (type == "p")
        out3d <- points3d
    if (!add)
        open3d()
    
    k <- dim(matrix)[1]
    if (is.null(size)) {
        if (type != "p")
            sz <- (cSize(matrix)/sqrt(k))*(1/80)
        else
            sz <- 10
    } else {
        sz <- size
    }
    if (size > 0)
        out3d(matrix,col=col1,radius=sz, size=sz)
    if(both) {
        if (size > 0)
            out3d(tarmatrix,col=col2,radius=sz, size=sz)
        if (lines) {
            linemesh <- list()
            linemesh$vb <- t(cbind(rbind(matrix,tarmatrix),1))
            linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
            class(linemesh) <- "mesh3d"
            wire3d(linemesh,lwd=1.5,col=lcol,lit=FALSE)
        }
    }
    x2 <- x1 <- x3 <- c(0:(ngrid-1))/ngrid;
    x0 <- as.matrix(expand.grid(x1,x2,x3))
    
    cent.mat <- apply(matrix,2,scale,scale=F)
    mean.mat <- apply(matrix,2,mean)
    
    if (ngrid > 1) {
        xrange <- diff(range(matrix[,1]))
        yrange <- diff(range(matrix[,2]))
        zrange <- diff(range(matrix[,3]))
        xrange1 <- diff(range(tarmatrix[,1]))
        yrange1 <- diff(range(tarmatrix[,2]))
        zrange1 <- diff(range(tarmatrix[,3]))
        maxi <- max(c(xrange,yrange,zrange,xrange1,yrange1,zrange1))
        maxi <- maxi+0.02*maxi
        x0 <- maxi*x0
        x0 <- apply(x0,2,scale,scale=FALSE)
        space <- eigen(crossprod(cent.mat))$vectors
        x0 <- t(t(x0%*%space)+mean.mat)
        x0 <- tps3d(x0,matrix,tarmatrix)
        
        ## create deformation cube
        outmesh <- list(vb = rbind(t(x0),1))
        class(outmesh) <- "mesh3d"
        
        yinit <- xinit <- zinit <- NULL
        if (2 %in% showaxis) {
            xinit0 <- xinit <- (c(1,2,2+ngrid,1+ngrid))
            for (i in 1:(ngrid-2))
                xinit <- cbind(xinit,(xinit0+i))
            
            xinit0 <- xinit
            for (i in 1:(ngrid-2))
                xinit <- cbind(xinit,xinit0+(i*ngrid))
            
            xinit0 <- xinit
            for (i in 1:(ngrid-1))
                xinit <- cbind(xinit,xinit0+(i*ngrid^2))
        }
        if (1 %in% showaxis) {
            yinit0 <- yinit <- c(ngrid,ngrid+ngrid^2, 2*ngrid+ngrid^2, 2*ngrid)
            for( i in 1:(ngrid-2))
                yinit <- cbind(yinit,yinit0+i*ngrid)
            yinit0 <- yinit
            for (i in 1:(ngrid-2))
                yinit <- cbind(yinit,yinit0+i*ngrid^2)
            yinit0 <- yinit
            for (i in 1:(ngrid-1))
                yinit <- cbind(yinit,yinit0-i)
        }
        if (3 %in% showaxis) {
            zinit0 <- zinit <- (c(2,1,1+ngrid^2,2+ngrid^2))
            for( i in 1:(ngrid-2))
                zinit <- cbind(zinit,(zinit0+i))
            
            zinit0 <- zinit
            for (i in 1:(ngrid-2))
                zinit <- cbind(zinit,zinit0+(i*ngrid^2))
            
            zinit0 <- zinit
            
            for (i in 1:(ngrid-1))
                zinit <- cbind(zinit,zinit0+(i*ngrid))
        }
        outmesh$ib <- cbind(xinit,yinit,zinit)
        wire3d(outmesh,lit=F)
    }
}


