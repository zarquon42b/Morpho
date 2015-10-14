#' calculate an affine transformation matrix
#'
#' calculate an affine transformation matrix
#' @param x fix landmarks
#' @param y moving landmarks
#' @param type set type of affine transformation: options are  "rigid", "similarity" (rigid + scale) and "affine",
#' @param reflection logical: if TRUE "rigid" and "similarity" allow reflections.
#' @param lambda numeric: regularisation parameter of the TPS.
#' @details
#' \code{x} and \code{y} can also be a pair of meshes with corresponding vertices.
#' @return returns a 4x4 (3x3 in 2D case)  transformation matrix or an object of class "tpsCoeff" in case of type="tps".
#' @note all lines containing NA, or NaN are ignored in computing the transformation.
#' @examples
#' data(boneData)
#' trafo <- computeTransform(boneLM[,,1],boneLM[,,2])
#' transLM <- applyTransform(boneLM[,,2],trafo)
#' @export
computeTransform <- function(x,y,type=c("rigid","similarity","affine","tps"),reflection=FALSE,lambda=1e-8) {
    if (inherits(x,"mesh3d"))
        x <- vert2points(x)
     if (inherits(y,"mesh3d"))
        y <- vert2points(y)
    type <- substr(type[1],1L,1L)
    ##check for missing entries
    xrows <- rowSums(x)
    yrows <- rowSums(y)
    xbad <- which(as.logical(is.na(xrows) + is.nan(xrows)))
    ybad <- which(as.logical(is.na(yrows) + is.nan(yrows)))
    bad <- unique(c(xbad,ybad))
    if (length(bad)) {
        message("some landmarks are missing and ignored for calculating the transform")
        x <- x[-bad,]
        y <- y[-bad,]
    }
    if (type %in% c("r","s")) {
        scale <- TRUE
        if (type == "r")
            scale <- FALSE
        trafo <- getTrafo4x4(rotonto(x,y,scale = scale,reflection=reflection))
    } else if (type=="a"){
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
    } else if (type == "t") {
        m <- ncol(y)
        Lall <- CreateL(y,lambda=lambda, output="Linv")
        Linv <- Lall$Linv
        m2 <- rbind(x,matrix(0,m+1,m))
        coeff <- as.matrix(Linv%*%m2)
        trafo <- list(refmat=y,tarmat=x,coeff=coeff,lambda=lambda)
        class(trafo) <- "tpsCoeff"
    } else {
        stop("Unknown transformation type")
    }
        
    return(trafo)
}
