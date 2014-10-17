#' apply affine transformation to data
#'
#' apply affine transformation to data
#' @param x matrix or mesh3d
#' @param trafo 4x4 transformation matrix (for mesh3d the matrix will be transformed to a 4x4 matrix)
#' @param inverse logical: if TRUE, the inverse of the transformation is applied
#' @return the transformed object
#' @examples
#' data(boneData)
#' rot <- rotonto(boneLM[,,1],boneLM[,,2])
#' trafo <- getTrafo4x4(rot)
#' boneLM2trafo <- applyTransform(boneLM[,,2],trafo)
#' @rdname applyTransform
#' @export
applyTransform <- function(x,trafo,inverse)UseMethod("applyTransform")

#' @rdname applyTransform
#' @export
applyTransform.matrix <- function(x,trafo,inverse=FALSE) {
    if (inverse)
        trafo <- solve(trafo)
    out <-homg2mat(trafo%*%mat2homg(x))
    return(out)
}

#' @rdname applyTransform
#' @export
applyTransform.mesh3d <- function(x,trafo,inverse=FALSE) {
    if (ncol(trafo) == 3)
        trafo <- mat3x3tomat4x4(trafo)
     if (inverse)
         trafo <- solve(trafo)
     x$vb <- trafo%*%x$vb
     if (!is.null(x$normals))
         x <- vcgUpdateNormals(x)
     return(x)
 }


mat3x3tomat4x4 <- function(x) {
    n <- ncol(x)
    x <- rbind(cbind(x,0),0);x[n+1,n+1] <-1
    return(x)
}
