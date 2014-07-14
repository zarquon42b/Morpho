#' rotate matrix of landmarks
#' 
#' rotate matrix of landmarks by using a rotation determined by two matrices.
#' 
#' A matrix is rotated by rotation parameters determined by two different
#' matrices. This is usefull, if the rotation parameters are to be estimated by
#' a subset of landmark coordinates.
#' 
#' @param X Matrix to be rotated
#' @param refmat reference matrix used to estimate rotation parameters
#' @param tarmat target matrix used to estimate rotation parameters
#' @param scale logical: requests scaling to minimize sums of squared distances
#' @param reflection logical: if TRUE, reflections are allowed.
#' @param weights vector of length k, containing weights for each landmark.
#' @param centerweight logical: if weights are defined and centerweigths=TRUE,
#' the matrix will be centered according to these weights instead of the
#' barycenter.
#' @param getTrafo logical: if TRUE, a 4x4 transformation matrix will also be returned.
#' @return if \code{getTrafo=FALSE} the transformed X will be returned,
#' else alist containing:
#' \item{Xrot}{the transformed matrix X}
#' \item{trafo}{a 4x4 transformation matrix}
#' @author Stefan Schlager
#' @seealso \code{\link{rotonto}},\code{\link{rotmesh.onto}}
#' 
#' @examples
#' 
#' 
#' data(nose)
#' shortnose.rot <-
#' rotonmat(shortnose.lm,shortnose.lm[1:9,],longnose.lm[1:9,])
#' 
#' ##view result
#' \dontrun{
#' deformGrid3d(shortnose.rot,shortnose.lm,ngrid=0)
#' }
#' @export
rotonmat <- function(X,refmat,tarmat,scale=TRUE,reflection=FALSE, weights=NULL, centerweight=FALSE,getTrafo=FALSE) {	
    ro <- rotonto(tarmat,refmat,scale=scale,signref=F,reflection=reflection, weights=weights, centerweight=centerweight)
    hmat <- getTrafo4x4(ro)
    Xrot <- homg2mat(hmat%*%mat2homg(X))
    if (!getTrafo)
        return(Xrot)
    else
        return(list(Xrot=Xrot,trafo=hmat))
}
