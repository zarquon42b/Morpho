#' Compute area enclosed within an irregular polygon
#'
#' Compute area enclosed within an irregular polygon - i.e. defined by curves
#' @param x k x 2 or k x 3 matrix containing ordered coordinates forming the boundary of the area. For 3D-cases, the area should be closed to a 2D surface (see details below). 
#' @details
#' For 3D coordinates, a PCA is computed and only the first two PCs are used to compute the area. This is a projection of the coordinates onto a 2D plane spanned by those PCs.
#' @note in case custom planes are preferred, the data can first be projected onto such a custom defined plane via \code{\link{points2plane}} first.
#' @return
#' returns a list containing
#' \item{area}{size of the enclosed area}
#'  \item{xpro2D}{projected coordinates of x in the 2D plane.}
#'  \item{poly}{object of class \code{sp} as defined by the \code{sp} package.}
#' \item{xpro3D}{For 3D-cases, this contains the projected coordinates of x rotated back into the original coordinate system}
#' @examples
#' require(shapes)
#' require(rgeos)
#' myarea <- computeArea(gorf.dat[c(1,6:8,2:5),,1])
#' myarea$area
#' plot(myarea$poly)
#'
#'
#' ## 3D example
#' data(boneData)
#' myarea3D <- computeArea(boneLM[c(4,2,3,7,5,6,8),,1])
#' plot(myarea3D$poly)
#' cent <- colMeans(myarea3D$xpro2D)
#' text(cent[1],cent[2],labels=paste0("Area=",round(myarea3D$area,digits=2)))
#' @export
computeArea <- function(x) {
    d <- ncol(x)
    k <- nrow(x)
    if (d == 3) {
    xpca <- Morpho::prcompfast(x)
    ## only use the two most important functions
    xpro <- xpca$x[,1:2]
   
    } else
        xpro <- x
    ### setup data to be evaluated by rgeos
    curveCat <- t(cbind(xpro,","))
    rownames(curveCat) <- NULL
    mycat <- c(curveCat,xpro[1,])
    mycat <- c("POLYGON((",mycat,"))")
    mycat1 <- paste0(mycat,collapse = " ")
   
    poly <- rgeos::readWKT(mycat1)
    area <- rgeos::gArea(poly)
    xpro3D <- NULL
    if (d == 3)
        xpro3D <- sweep(xpro %*%t(xpca$rotation[,1:2]),2,-xpca$center[1:3])
    return(list(area=area,xpro2D=xpro,poly=poly,xpro3D=xpro3D))
    
}
