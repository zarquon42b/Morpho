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
#' @param \dots Additional parameters which will be passed to the methods.
#' @return returns an invisible array containing the shapes associated with the Principal components selected.
#' @seealso \code{\link{procSym}}
#' @examples
#' 
#' \dontrun{
#' data(nose)
#' #make a tiny sample
#' nosearr <- bindArr(longnose.lm, shortnose.lm, along=3)
#' proc <- procSym(nosearr)
#' pcaplot3d(proc,pcshow=1,mag=-3)#only one PC available
#' }
#' @rdname pcaplot3d
#' @export
pcaplot3d <- function (x,...) UseMethod("pcaplot3d")

#' @rdname pcaplot3d
#'
#' @export
pcaplot3d.symproc <- function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,sym=TRUE,...) 
{   
  refshape <- x$mshape
  if (sym) {
      PCs <- x$PCsym
      Scores <- x$PCscore_sym
  } else {
      PCs <- x$PCasym
      Scores <- x$PCscore_asym
  }
  
  .pcaplot3d(refshape, PCs, Scores, pcshow=pcshow, mag=mag,color=color, lwd=lwd)
}
#' @rdname pcaplot3d
#'
#' @export
pcaplot3d.nosymproc <- function(x,pcshow=c(1,2,3),mag=3,color=4,lwd=1,...)
{   
    refshape <- x$mshape
    PCs <- x$PCs
    Scores <- x$PCscores
    .pcaplot3d(refshape, PCs, Scores, pcshow=pcshow, mag=mag,color=color, lwd=lwd)
}

.pcaplot3d <- function(refshape,PCs, Scores, pcshow=c(1,2,3), mag=3,color=4,lwd=1,...) {
    A <- refshape
    k <- dim(A)[1]
    m <- dim(A)[2]
    if (is.vector(PCs))
        PCs <- matrix(PCs,length(PCs),1)
    
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
    spheres3d(refshape,  col = color,radius=sz)
    invisible(outarr)
}
