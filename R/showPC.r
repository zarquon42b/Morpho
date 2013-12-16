#' convert PCs to landmark configuration
#' 
#' convert PC-scores to landmark coordinates
#' 
#' Rotates and translates PC-scores derived from shape data back into
#' configuration space.
#' 
#' @param scores vector of PC-scores
#' @param PC Principal components (eigenvectors of the covariance matrix)
#' associated with 'scores'.
#' @param mshape matrix containing the meanshape's landmarks (used to center
#' the data by the PCA)
#' @return returns matrix containing landmarks
#' @author Stefan Schlager
#' @seealso \code{\link{prcomp}}, \code{\link{procSym}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' library(shapes)
#' ## generate landmarks using
#' ##the first PC-score of the first specimen
#' 
#' proc <- procSym(gorf.dat)
#' lm <- showPC(proc$PCscores[1,1],proc$PCs[,1],proc$mshape)
#' plot(lm,asp=1)
#' 
#' ##now the first 3 scores
#' lm2 <- showPC(proc$PCscores[1,1:3],proc$PCs[,1:3],proc$mshape)
#' points(lm2,col=2)
#' 
#' @export
showPC <- function(scores,PC,mshape)
  {
    dims <- dim(mshape)
    if (length(scores) > 1)                                    #pred1 <- t(coeff)%*%mod
      predPC <- PC%*%scores
    else
      predPC <- PC*scores
    modell <- mshape+matrix(predPC,dims[1],dims[2])
    return(modell)
}
