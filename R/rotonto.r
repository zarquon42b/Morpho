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
rotonto <- function(x,y,scale=FALSE,signref=TRUE,reflection=TRUE,weights=NULL,centerweight=FALSE) {
    reflect=0
    m <- dim(x)[2]
    if (!is.null(weights))
        weights <- weights/sum(weights)
    
    X <- apply(x,2,scale,scale=F)
    Y <- apply(y,2,scale,scale=F)
    if (centerweight && !is.null(weights)) {
        xcent <- apply(X*weights,2,sum)
        ycent <- apply(Y*weights,2,sum)
        X <- scale(X,scale=F,center=xcent)
        Y <- scale(Y,scale=F,center=ycent)
    }
    if (!is.null(weights)) {
        Dn <- diag(weights)
        X1 <- Dn%*%X
        Y1 <- Dn%*%Y
        XY <- crossprod(X1,Y1)
    } else {
        XY <- crossprod(X,Y)
    }
    sv1 <- svd(XY)
    gamm <- tcrossprod(sv1$v,sv1$u)
    
    if(sign(det(gamm))<1)
        {	reflect <- 1
		if (signref && reflection)
                    cat("reflection involved\n")
                
                if (!reflection) {
                    u <- sv1$u
                    v <- sv1$v
                    chk1 <- Re(prod(eigen(v)$values))
                    chk2 <- Re(prod(eigen(u)$values))
                    if ((chk1 < 0) && (chk2 > 0)) {
                        v[, dim(v)[2]] <- v[, dim(v)[2]] * (-1)
                        gamm <- v %*% t(u)
                    }
                    if ((chk2 < 0) && (chk1 > 0)) {
                        u[, dim(u)[2]] <- u[, dim(u)[2]] * (-1)
                        gamm <- v %*% t(u)
                    }
                }
            }
    trans <- x[1,]-X[1,]
    transy <- y[1,]-Y[1,]
    del <- sv1$d
    
    ctrace <- function(MAT) sum(diag(crossprod(MAT)))
    if (scale) {
        if (!is.null(weights))
            bet <- sum(del)/ctrace(Y1)
        else
            bet <- sum(del)/ctrace(Y)
        yrot <- bet*Y%*%gamm
    } else {
        bet <- 1
        yrot <- Y%*%gamm
    }
    Y <- yrot  	
    yrot <- t(t(yrot)+trans)
    out <- list(yrot=yrot,Y=Y,X=X,trans=trans,transy=transy,gamm=gamm,bet=bet,reflect=reflect)
    class(out) <- "rotonto"
    return(out)
}

#' @rdname rotonto
#' @export
rotreverse <- function(mat,rot)UseMethod("rotreverse")

#' @rdname rotonto
#' @export
rotreverse.matrix <- function(mat,rot){
    hmat <- solve(getTrafo4x4(rot))
    out <-homg2mat(hmat%*%mat2homg(mat))
    return(out)
}

#' @rdname rotonto
#' @export
rotreverse.mesh3d <- function(mat,rot)
    {
        x <- rotreverse(vert2points(mat),rot)
        mat$vb[1:3,] <- t(x)
         if (!is.null(mat$normals))
        mat <- updateNormals(mat)

        return(mat)
    }

#' get 4x4 Transformation matrix
#'
#' get 4x4 Transformation matrix
#' @param object of class "rotonto"
#'
#' @return returns a 4x4 transformation matrix
#' @examples
#' data(boneData)
#' rot <- rotonto(boneLM[,,1],boneLM[,,2])
#' trafo <- getTrafo4x4(rot)
#' @rdname getTrafo4x4
#' @export
getTrafo4x4 <- function(x)UseMethod("getTrafo4x4")

#' @rdname getTrafo4x4
#' @export
getTrafo4x4.rotonto <- function(x) {
    m <- ncol(x$gamm)
    hgamm <- rbind(cbind(x$gamm,0),0);hgamm[m+1,m+1] <- 1
    htrans <- diag(m+1);htrans[1:m,m+1] <- c(-x$transy)
    htrans2 <- diag(m+1);htrans2[1:m,m+1] <- c(x$trans)
    scale <- diag(m+1);diag(scale)[1:m] <- x$bet
    hall <- htrans2%*%scale%*%t(hgamm)%*%htrans
    return(hall)
}

mat2homg <- function(x) {
    x <- rbind(t(x),1)
    return(x)
}

homg2mat <- function(x) {
    x <- t(x[1:3,])
    return(x)
}
#' apply affine transformation to data
#'
#' apply affine transformation to data
#' @param x matrix or mesh3d
#' @param trafo 4x4 transformation matrix
#' @param inverse logical: if TRUE, the inverse of the transformation is applied
#' @return the transformed object
#' @examples
#' data(boneData)
#' rot <- rotonto(boneLM[,,1],boneLM[,,2])
#' trafo <- getTrafo4x4(rot)
#' boneLM2trafo <- applyTransformation(boneLM[,,2],trafo)
#' @rdname applyTransformation
#' @export
applyTransformation <- function(x,trafo,inverse)UseMethod("applyTransformation")

#' @rdname applyTransformation
#' @export
applyTransformation.matrix <- function(x,trafo,inverse=FALSE) {
    if (inverse)
        trafo <- solve(trafo)
    out <-homg2mat(trafo%*%mat2homg(x))
    return(out)
}

#' @rdname applyTransformation
#' @export
applyTransformation.mesh3d <- function(x,trafo,inverse=FALSE) {
     if (inverse)
         trafo <- solve(trafo)
     x$vb <- trafo%*%x$vb
     if (!is.null(x$normals))
         x <- updateNormals(x)
     return(x)
 }
