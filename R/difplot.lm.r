.difplotLM <- function(refshape,targetshape,color=4,lwd=1,lcol=2,rgl.new=TRUE, text=TRUE)
{   
    if (rgl.new == TRUE)
        open3d()
    
    A <- refshape
    k <- dim(A)[1]
    m <- dim(A)[2]
    sds <- 0
     lim <- max(abs(refshape))
    sz <- (cSize(refshape)/sqrt(k))*(1/80)
    spheres3d(refshape,  col = color,radius=sz)
    linemesh <- list()
    linemesh$vb <- t(cbind(rbind(refshape,targetshape),1))
    linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
    class(linemesh) <- "mesh3d"
    wire3d(linemesh,col=lcol,lwd=lwd,lit=F)
    if (text)
        text3d(refshape,texts=paste("",c(1:k),sep=""),cex=1,col=lcol,adj=1.2) 
}
