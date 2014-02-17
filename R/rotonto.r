#' rotates, translates and scales one matrix onto an other using Procrustes
#' fitting
#' 
#' rotate a matrix onto an other without loosing information about the location
#' of the targetmatrix and reverse this transformations using rotreverse
#' 
#' 
#' @title rotates, translates and scales one matrix onto an other using Procrustes
#' fitting
#' @param x k x m matrix to be rotated onto (targetmatrix)
#' @param y k x m matrix which will be rotated (reference matrix)
#' @param scale logical: scale matrix to minimize sums of squares
#' @param signref logical: report if reflections were involved in the rotation
#' @param mat matrix on which the reverse transformations have to be applied
#' @param rot an object resulting from the former application of rotonto
#' @param reflection allow reflections.
#' @param weights vector of length k, containing weights for each landmark.
#' @param centerweight logical: if weights are defined and centerweigths=TRUE,
#' the matrix will be centered according to these weights instead of the
#' barycenter.
#' @return
#' \item{yrot }{rotated and translated matrix}
#' \item{Y }{centred and rotated reference matrix}
#' \item{X }{centred target matrix}
#' \item{trans }{vector between original position of target and centered
#' reference (during rotation process)}
#' \item{transy }{vector between original position of reference and
#' centered reference (during rotation process)}
#' \item{gamm }{rotation matrix}
#' \item{bet }{scaling factor applied}
#' \item{reflect }{if \code{reflect = 1}, reflections are involved in the
#' superimposition. Else, reflect = 0}
#' @author Stefan Schlager
#' @seealso \code{\link{rotmesh.onto}}
#' @references Lissitz, R. W., Sch6nemann, P. H., & Lingoes, J. C. (1976). A
#' solution to the weighted Procrustes problem in which the transformation is
#' in agreement with the loss function. Psychometrika, 41,547-550.
#' 
#' @examples
#' 
#' library(shapes)
#' lims <- c(min(gorf.dat[,,1:2]),max(gorf.dat[,,1:2]))
#' rot <- rotonto(gorf.dat[,,1],gorf.dat[,,2]) ### rotate the second onto the first config
#' plot(rot$yrot,pch=19,xlim=lims,ylim=lims) ## view result
#' points(gorf.dat [,,2],pch=19,col=2) ## view original config
#' rev1 <- rotreverse(rot$yrot,rot)
#' points(rev1,cex=2) ### show inversion by larger circles around original configuration
#' 
#' 
#' @export
rotonto <- function(x,y,scale=FALSE,signref=TRUE,reflection=TRUE,weights=NULL,centerweight=FALSE)
{ 	reflect=0
  	m <- dim(x)[2]
        if (!is.null(weights))
          weights <- weights/sum(weights)
       
        X <- apply(x,2,scale,scale=F)
  	Y <- apply(y,2,scale,scale=F)
         if (centerweight && !is.null(weights))
          {
            
            xcent <- apply(X*weights,2,sum)
            ycent <- apply(Y*weights,2,sum)
            X <- scale(X,scale=F,center=xcent)
            Y <- scale(Y,scale=F,center=ycent)
          }
        if (!is.null(weights))
          {
            Dn <- diag(weights)
            X1 <- Dn%*%X
            Y1 <- Dn%*%Y
            XY <- crossprod(X1,Y1)
          }
        else
          XY <- crossprod(X,Y)

        sv1 <- svd(XY)
        
	#dd <- diag(sign(sv1$d))
	gamm <- tcrossprod(sv1$v,sv1$u)
	#gamm <- (sv1$v)%*%gamm
  	
  	if(sign(det(gamm))<1)
	  {	reflect <- 1
		if (signref && reflection)
                  {cat("reflection involved\n")
                 }
                if (!reflection)
                  {
                    u <- sv1$u
                    v <- sv1$v
                    chk1 <- Re(prod(eigen(v)$values))
                    chk2 <- Re(prod(eigen(u)$values))
                    if ((chk1 < 0) && (chk2 > 0))
                      {
                        v[, dim(v)[2]] <- v[, dim(v)[2]] * (-1)
                        gamm <- v %*% t(u)
                      }
                    if ((chk2 < 0) && (chk1 > 0))
                      {
                        u[, dim(u)[2]] <- u[, dim(u)[2]] * (-1)
                        gamm <- v %*% t(u)
                      }
                  }
              }
                
                  
	  
  	trans <- x[1,]-X[1,]
  	transy <- y[1,]-Y[1,]
  	#yrot <- Y%*%gamm
  	#sig <- sign(det(XY))
  	del <- sv1$d
   #del[m] <- sig*abs(del[m])
  	ctrace <- function(MAT) sum(diag(crossprod(MAT)))
  
  
  	if (scale == TRUE)
    		{
                  if (!is.null(weights))
                    bet <- sum(del)/ctrace(Y1)
                  else
                    bet <- sum(del)/ctrace(Y)
                  yrot <- bet*Y%*%gamm
		}
  	else
    		{
                  bet <- 1
                  yrot <- Y%*%gamm
		}
	Y <- yrot  	
	yrot <- t(apply(yrot,1,function(x){x+trans}))
  
  	return(list(yrot=yrot,Y=Y,X=X,trans=trans,transy=transy,gamm=gamm,bet=bet,reflect=reflect))
}
#' @rdname rotonto
#' @export
rotreverse <- function(mat,rot)
  {
    transfun <- function(x,trans)
      {
        x <- x+trans
      }
    
    out <- t(apply(mat,1,transfun,-rot$trans))
    out <- t(apply((out%*%t(rot$gamm))*1/rot$bet,1,transfun,rot$transy))
    return(out)
  }
