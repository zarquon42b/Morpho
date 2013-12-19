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
#' @return returns rotated X.
#' @author Stefan Schlager
#' @seealso \code{\link{rotonto}},\code{\link{rotmesh.onto}}
#' @keywords ~kwd1 ~kwd2
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
rotonmat <- function(X,refmat,tarmat,scale=TRUE,reflection=FALSE, weights=NULL, centerweight=FALSE)
{	
	ro <- rotonto(tarmat,refmat,scale=scale,signref=F,reflection=reflection, weights=weights, centerweight=centerweight)
	
	Xrot <- t(apply(X,1,function(x){x-ro$transy}))%*%ro$gamm
	if (scale)
		{sf <- cSize(refmat)/cSize(ro$yrot)
		Xrot <- Xrot/sf
		}
	Xrot <- t(apply(Xrot,1,function(x){x+ro$trans}))
	return(Xrot)
}
