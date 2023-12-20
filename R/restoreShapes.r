#' restore shapes from PC-Scores or similar projections
#'
#' restore shapes from PC-Scores or similar projections
#' 
#' Rotates and translates PC-scores (or similar) derived from shape data back into
#' configuration space.
#' 
#' @param scores vector of PC-scores, or matrix with rows containing PC-scores
#' @param PC Principal components (eigenvectors of the covariance matrix)
#' associated with 'scores'.
#' @param mshape matrix containing the meanshape's landmarks (used to center
#' the data by the PCA)
#' @param sizeshape logical: if TRUE, it is assumed that the data is the output of \code{procSym} run with \code{sizeshape=TRUE}.
#' @param origsize logical: if \code{sizeshape = TRUE}, this will apply the scaling to the original size from the corresponding entry from the PC basis matrix.
#' @param meanlogCS numeric: provide the average log Centroid Size of the original sample (see examples below). Only needed if \code{sizeshape = TRUE} and \code{origsize = TRUE}
#' @return returns matrix or array containing landmarks
#' @author Stefan Schlager
#' @seealso \code{\link{prcomp}}, \code{\link{procSym}}
#' 
#' @examples
#' 
#' if (require(shapes)) {
#' ## generate landmarks using
#' ##the first PC-score of the first specimen
#' 
#' proc <- procSym(gorf.dat)
#' lm <- restoreShapes(proc$PCscores[1,1],proc$PCs[,1],proc$mshape)
#' plot(lm,asp=1)
#' 
#' ##now the first 3 scores
#' lm2 <- restoreShapes(proc$PCscores[1,1:3],proc$PCs[,1:3],proc$mshape)
#' points(lm2,col=2)
#'
#' ## Now restore some sizeshape data
#' procSize <- procSym(gorf.dat,sizeshape=TRUE)
#' est1 <- restoreShapes(range(procSize$PCscores[,1]),procSize$PCs[,1],procSize$mshape,
#'                       sizeshape=TRUE,origsize=TRUE,meanlogCS=procSize$meanlogCS)
#' }
#' 
#' @seealso \code{\link{getPCscores}}
#' @export
restoreShapes <- function(scores,PC,mshape,sizeshape=FALSE,origsize=FALSE,meanlogCS)
  {
    dims <- dim(mshape)
    PC <- as.matrix(PC)
    
    if (!is.matrix(scores) && ncol(PC) == 1)
        if (length(scores) > 1)
            scores <- as.matrix(scores)
    if (!is.matrix(scores)){
        if (length(scores) != ncol(PC))
            stop("scores must be of the same length as ncol(PC)")
        predPC <- PC%*%scores
        if (!sizeshape)
            modell <- mshape+matrix(predPC,dims[1],dims[2])
        else {
            modell <- mshape+matrix(predPC[-1],dims[1],dims[2])
            if (origsize) {
                
                if (missing(meanlogCS))
                    stop("please provide mean log centroid size")
            
                modell <- modell*(exp(predPC[1]+meanlogCS))
            }
        }
        return(modell)
    } else {
          n <- nrow(scores)
          outarr <- array(0,dim=c(dims,n))
          for (i in 1:n) {
              outarr[,,i] <- restoreShapes(scores[i,],PC,mshape,sizeshape=sizeshape,origsize=origsize,meanlogCS=meanlogCS)
          }
          if (!is.null(rownames(scores)))
              dimnames(outarr)[[3]] <- rownames(scores)
          return(outarr)
    }    
}

#' Obtain PC-scores for new landmark data
#'
#' Obtain PC-scores for new landmark data
#' @param x landmarks aligned (e.g. using \code{\link{align2procSym}} to the meanshape of data the PCs are derived from.
#' @param PC Principal components (eigenvectors of the covariance matrix)
#' @param mshape matrix containing the meanshape's landmarks (used to center the data)
#' @return returns a matrix containing the PC scores
#' @examples
#' \dontrun{
#' data(boneData)
#' proc <- procSym(boneLM[,,-c(1:2)])
#' newdata <- boneLM[,,c(1:2)]
#' newdataAlign <- align2procSym(proc,newdata)
#' scores <- getPCscores(newdataAlign,proc$PCs,proc$mshape)
#' }
#' @seealso \code{\link{restoreShapes}}
#' @export
getPCscores <- function(x, PC, mshape) {
    if (is.matrix(x))
        x <- array(x,dim=(c(dim(x),1)))
    x <- sweep(x,1:2,mshape)
    x <- vecx(x)
    scores <- x%*%(PC)#%*%t(x)
    return(scores)
}

#' restore original data from PCA
#'
#' restore original data from PCA by reverting rotation and centering
#' @param scores matrix containing the PC-scores
#' @param rotation matrix containing the PCs
#' @param center vector containing the center
#'
#' @examples
#' myirispca <- prcomp(iris[,1:4])
#' myirisRecovered <- restoreFromPCA(myirispca$x,myirispca$rotation,myirispca$center)
#' all.equal(myirisRecovered,as.matrix(iris[,1:4]))
#' @export
restoreFromPCA <- function(scores,rotation,center) {
    predPC <- t(as.matrix(rotation) %*% t(scores))
    predPC <- sweep(predPC,2,-center)
}
