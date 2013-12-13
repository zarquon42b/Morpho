#' @rdname render
#' @export
export <- function(x,...)UseMethod("export")

#' @rdname render
#' @method export meshDist
#' @export
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
