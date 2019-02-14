#' visualise differences between two superimposed sets of 2D landmarks
#' 
#' visualise differences between two superimposed sets of 2D landmarks by
#' deforming a square grid based on a thin-plate spline interpolation
#' 
#' 
#' @param matrix reference matrix containing 2D landmark coordinates or mesh of class "mesh3d"
#' @param tarmatrix target matrix containing 2D landmark coordinates or mesh of class "mesh3d"
#' @param ngrid number of grid lines to be plotted; ngrid=0 suppresses grid
#' creation.
#' @param lwd width of lines connecting landmarks.
#' @param show integer (vector): if c(1:2) both configs will be plotted, show = 1 only plots the reference and show = 2 the target.
#' plotted. Options are combinations of 1,2 and 3.
#' @param lines logical: if TRUE, lines between landmarks will be plotted.
#' @param lcol color of lines
#' @param col1 color of "matrix"
#' @param col2 color of "tarmat"
#' @param pcaxis logical: align grid by shape's principal axes.
#' @param add logical: if TRUE, output will be drawn on existing plot.
#' @param wireframe list/vector containing row indices to be plotted as wireframe (see \code{\link{lineplot}}.)
#' @param margin margin around the bounding box to draw the grid
#' @param gridcol color of the grid
#' @param cex1 control size of points belonging to \code{matrix}
#' @param cex2 control size of points belonging to \code{tarmatrix}
#' @param ... additional parameters passed to plot
#' @author Stefan Schlager
#' @seealso \code{\link{tps3d}}
#' 
#' @examples
#' if (require(shapes)) {
#' proc <- procSym(gorf.dat)
#' deformGrid2d(proc$mshape,proc$rotated[,,1],ngrid=5,pch=19)
#' }
#' 
#' @export

deformGrid2d <- function(matrix,tarmatrix,ngrid=0,lwd=1,show=c(1:2),lines=TRUE,lcol=1,col1=2,col2=3,pcaxis=FALSE,add=FALSE,wireframe=NULL,margin=0.2,gridcol="grey",cex1=1,cex2=1,...)
{
    k <- dim(matrix)[1]
    x0 <- NULL
    if (ngrid > 1) {

        x2 <- x1 <- c(0:(ngrid-1))/ngrid;
        x0 <- as.matrix(expand.grid(x1,x2))
        
       
        xrange <- diff(range(matrix[,1]))
        yrange <- diff(range(matrix[,2]))
        
        xrange1 <- diff(range(tarmatrix[,1]))
        yrange1 <- diff(range(tarmatrix[,2]))
        
        ## cent.mat <- scale(matrix,scale=FALSE)
        ## mean.mat <- colMeans((matrix+tarmatrix)/2)
        mean.mat <- c((min(matrix[,1])+max(matrix[,1]))/2,(min(matrix[,2])+max(matrix[,2]))/2)
        cent.mat <- scale(matrix,scale=FALSE,center=mean.mat)
        maxi <- max(c(xrange,yrange,xrange1,yrange1))
        maxi <- (1+margin)*maxi
        x0 <- maxi*x0
        x0 <- scale(x0, scale=FALSE)
        x0[,2] <- x0[,2]
        if (pcaxis)
            space <- eigen(crossprod(cent.mat))$vectors
        else
            space <- diag(2)
        
        x0 <- (x0%*%space)
        x00 <- x0 <- scale(x0,center=-mean.mat,scale=F)
        x0 <- tps3d(x0,matrix,tarmatrix,threads=1)
        
    }
    lims <- apply(rbind(matrix,tarmatrix,x0),2,range)
    if (1 %in% show) {
        if (add)
            points(matrix,col=col1,cex=cex1,...)
        else 
            plot(matrix,col=col1,xlim=lims[,1],ylim = lims[,2],asp=1,xlab = "",ylab = "", axes = F,cex=cex1,...)
        if (!is.null(wireframe))
            lineplot(matrix,wireframe,col=col1,lwd=lwd)
        
    }
    if(2 %in% show) {
        if (1 %in% show || add)
            points(tarmatrix,col=col2,,cex=cex2,...)
        else 
            plot(tarmatrix,col=col2,xlim=lims[,1],ylim = lims[,2],asp=1,xlab = "",ylab = "", axes = F,,cex=cex2,...)
        if (!is.null(wireframe))
            lineplot(tarmatrix,wireframe,col=col2,lwd=lwd)
    }
    if (lines) {
        linemesh <- list()
        linemesh$vb <- rbind(matrix,tarmatrix)
        linemesh$it <- cbind(1:k,(1:k)+k)
        for (i in 1:nrow(linemesh$it))
            lines(linemesh$vb[linemesh$it[i,],],lwd=lwd,col=3,lty=2)
    }
    
    if (ngrid > 1) {
        myrange <- 0:(ngrid-1)
        for (i in 0:(ngrid-1)) {
            lines(x0[(1:ngrid)+(i*ngrid),],col=gridcol)
            lines(x0[(myrange*ngrid)+i+1,],col=gridcol)
        }
    }
}

