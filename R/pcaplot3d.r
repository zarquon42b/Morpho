#' visualization of shape change
#' 
#' visualization of the shape changes explained by Principal components
#' 
#' @title visualization of shape variation
#' @param x a object derived from the function procSym calculated on 3D
#' coordinates.
#' @param pcshow a vector containing the PCscores to be visualized.
#' @param mag a vector or an integer containing which standard deviation of
#' which PC has to be visualized.
#' @param color color of the 3d points/spheres.
#' @param lwd width of the lines representing the shape change.
#' @param sym logical: if TRUE the symmetric component of shape is displayed.
#' Otherwise the asymmetric one.
#' @param legend logical: if TRUE a legend explaining the color coding of the PCs is plotted.
#' @param type character: for \code{type="spheres"}, the landmarks will be rendered using rgl's \code{spheres3d} function and for \code{type="points"} by  \code{points3d} respectivly.
#' @param \dots Additional parameters which will be passed to the methods.
#' @return returns an invisible array containing the shapes associated with the Principal components selected.
#' @seealso \code{\link{procSym}}
#' @examples
#' 
#' \dontrun{
#' data(boneData)
#' proc <- procSym(boneLM)
#' pcaplot3d(proc,pcshow=1:3,mag=-3)#only one PC available
#' }
#' @rdname pcaplot3d
#' @export
pcaplot3d <- function (x,...) UseMethod("pcaplot3d")

#' @rdname pcaplot3d
#'
#' @export
pcaplot3d.symproc <- function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,sym=TRUE,legend=TRUE,type=c("spheres","points"),...) 
{   
  refshape <- x$mshape
  if (sym) {
      PCs <- x$PCsym
      Scores <- x$PCscore_sym
  } else {
      PCs <- x$PCasym
      Scores <- x$PCscore_asym
  }
  
  .pcaplot3d(refshape, PCs, Scores, pcshow=pcshow, mag=mag,color=color, lwd=lwd,legend=legend,type=type)
}
#' @rdname pcaplot3d
#'
#' @export
pcaplot3d.nosymproc <- function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,legend=TRUE,type=c("spheres","points"),...)
{   
    refshape <- x$mshape
    PCs <- x$PCs
    Scores <- x$PCscores
    .pcaplot3d(refshape, PCs, Scores, pcshow=pcshow, mag=mag,color=color, lwd=lwd, legend=legend, type=type)
}

.pcaplot3d <- function(refshape,PCs, Scores, pcshow=c(1,2,3), mag=3,color=4,lwd=1,legend=TRUE,type=c("spheres","points"),...) {
    A <- refshape
    k <- dim(A)[1]
    m <- dim(A)[2]
    if (is.vector(PCs))
        PCs <- matrix(PCs,length(PCs),1)
    type <- match.arg(type,c("spheres","points"))
    npc <- dim(PCs)[2]
    lpc <- length(pcshow)
    rainb <- rainbow(lpc)
    sds <- 0
    if (length(mag)==1)
        mag <- c(rep(mag,lpc))
        
    for (i in 1:npc)
        sds[i] <- sd(Scores[,i])
    sz <- cSize(refshape)/sqrt(k)*(1/80)
    outarr <- array(NA, dim=c(dim(refshape),length(pcshow)))
    for (i in 1:length(pcshow)) {
        pc <- refshape+matrix(PCs[,pcshow[i]]*mag[i]*sds[pcshow[i]],k,3)
        outarr[,,i] <- pc
        linemesh <- list()
        linemesh$vb <- t(cbind(rbind(refshape,pc),1))
        linemesh$it <- t(cbind(1:k,1:k,(1:k)+k))
        class(linemesh) <- "mesh3d"
        wire3d(linemesh,lwd=lwd,lit=F,col=rainb[i])
    }
    if (legend) {
        plot(0,0, xlab="", ylab="", axes =F, cex=0,xlim=c(-1,1), ylim=c(-1,1))
        legend(-1,1, pch=20, cex=2, col=rainb, legend=paste0("PC",pcshow))
    }
    if (type == "spheres")
        spheres3d(refshape,  col = color,radius=sz)
    else
        points3d(refshape,col=color)
    invisible(outarr)
}
