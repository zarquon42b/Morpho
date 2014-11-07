#' calculate an affine transformation matrix
#'
#' calculate an affine transformation matrix
#' @param x fix landmarks
#' @param y moving landmarks
#' @param type set type of affine transformation: options are  "rigid", "similarity" (rigid + scale) and "affine",
#' @param reflection logical: if TRUE "rigid" and "similarity" allow reflections.
#' @return returns a 4x4 (3x3 in 2D case)  transformation matrix
#' 
#' @examples
#' data(boneData)
#' trafo <- computeTransform(boneLM[,,1],boneLM[,,2])
#' transLM <- applyTransform(boneLM[,,2],trafo)
#' @export
computeTransform <- function(x,y,type=c("rigid","similarity","affine"),reflection=FALSE) {
    type <- substr(type[1],1L,1L)
    if (type %in% c("r","s")) {
        scale <- TRUE
        if (type == "r")
            scale <- FALSE
        trafo <- getTrafo4x4(rotonto(x,y,scale = scale,reflection=reflection))
    } else {
        k <- nrow(x)
        m <- ncol(x)
        xp <- as.vector(t(x))
        yh <- cbind(y,1)
        M <- matrix(0,k*m,m*(m+1))
        M[(1:k)*m-(m-1),1:(m+1)] <- yh
        M[(1:k)*m-(m-2),(m+2):(2*(m+1))] <- yh
        if (m == 3)    
            M[(1:k)*3,(m+6):(m+9)] <- yh
        projS <- armaGinv(M) %*%xp
        trafo <- matrix(projS,m,m+1,byrow = T)
        trafo <- rbind(trafo,0)
        trafo[m+1,m+1] <- 1
    }
    return(trafo)
}
