#' Get a point along a line with a given distance from the start of the line
#'
#' Get a point along a line with a given distance from the start of the line
#'
#' @param mat matrix with rows containing sequential coordinates
#' @param dist numeric: distance from the first point on the line.
#' @param reverse logical: if TRUE start from the end of the line
#'
#' @return returns a vector containing the resulting coordinate
#' 
#' @export
getPointAlongOutline <- function(mat,dist=15,reverse=FALSE) {
    mat <- as.matrix(mat)
    if (reverse)
        mat <- mat[nrow(mat):1,]
    dists <- sapply(2:nrow(mat),function(x) x <- sqrt(sum((mat[x,]-mat[x-1,])^2)))
    cumdist <- cumsum(dists)
    myind <- which.max(cumdist[which(cumdist < dist)])
    if (length(myind)) {
        mydist <- cumdist[myind]
        final <- mat[myind+2,]-mat[myind+1,]
        final <- final/sqrt(sum(final^2))
        final.length <- dist-mydist
        final.point <- mat[myind+1,]+final.length*final
    } else {
        final <- mat[2,]-mat[1,]
        final <- final/sqrt(sum(final^2))
        final.point <- mat[1,]+final*dist
    }
    return(final.point)
}
geoDist <- function(mat) {
    dists <- sapply(2:nrow(mat),function(x) x <- sqrt(sum((mat[x,]-mat[x-1,])^2)))
    return(sum(dists))
}

#' Resample a curve equidistantly
#'
#' Resample a curve equidistantly (optionally with smoothing)
#' @param x matrix containing coordinates
#' @param n number of resulting points on the resampled curve
#' @param smooth logical: if TRUE, the resulting curve will be smoothed by using bezier curves.
#' @param smoothn integer: define the refinement of the bezier curve. The higher this value, the closer the final curve will be to the original.
#' @param open logical: define whether it is a closed curve or not.
#' @return returns a matrix containing the resampled curve
#' @examples
#' data(nose)
#' x <- shortnose.lm[c(304:323),]
#' xsample <- resampleCurve(x,n=50)
#' @export
#' @importFrom bezier bezier
resampleCurve <- function(x,n,smooth=FALSE,smoothn=n,open=TRUE) {
  
    if (!open) {
        x <- rbind(x,x[1,])
        n <- n+1
    }
    gd <- geoDist(x)
    dists <- seq(from=0,to=gd,length.out = n)
    

    if (smooth) {
        nan <- TRUE
        while(nan) {
            dists1 <- seq(from=0,to=gd,length.out = smoothn)
            out <- t(sapply(dists1,function(y) y <- t(getPointAlongOutline(x,dist=y))))
            t <- seq(0, 1, length=smoothn)
            out1 <- bezier(t,out)
            nan <- sum(is.nan(out1))
            if (nan) {
                smoothn <- floor(smoothn*.8)
                warning(paste0("bezier curve fit not successful: smoothn is set to ",smoothn))
            }
        }
        out <- resampleCurve(out1,n,smooth=FALSE)
    } else {
        out <- t(sapply(dists,function(y) y <- t(getPointAlongOutline(x,dist=y))))
    }
    if (!open)
        out <- out[-nrow(out),]
    
    return(out)
}

