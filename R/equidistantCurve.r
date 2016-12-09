#' make a curve equidistant (optionally up/downsampling)
#'
#' make a curve equidistant (optionally up/downsampling)
#' @param x k x m matrix containing the 2D or 3D coordinates
#' @param n integer: number of coordinates to sample. If NULL, the existing curve will be made equidistant.
#' @param open logical: specifies whether the curve is open or closed.
#' @param subsample integer: number of subsamples to draw from curve for interpolation. For curves with < 1000 points, no subsampling is required.
#' @param increment integer: if > 1, the curve is estimated iteratively by incrementing the original points by this factor. The closer this value to 1, the smoother the line but possibly farther away from the control points.
#' @param smoothit integer: smoothing iterations after each step
#' @param mesh specify mesh to project point to
#' @param iterations integer: how many iterations to run equidistancing.
#' @details
#' Equidistancy is reached by iteratively deforming (using TPS) a straight line with equidistantly placed points to the target using control points with the same spacing as the actual curve. To avoid singularity, the straight line containes a small amount of noise, which can (optionally) be accounted for by smoothing the line by its neighbours.
#' @note if n >> number of original points, the resulting curves can show unwanted distortions.
#' @return
#' matrix containing equidistantly placed points
#' @examples
#' data(nose)
#' x <- shortnose.lm[c(304:323),]
#' xsample <- equidistantCurve(x,n=50,iterations=10,increment=2)
#' \dontrun{
#' require(rgl)
#' points3d(xsample,size=5)
#' spheres3d(x,col=2,radius=0.3,alpha=0.5)
#' }
#' @export
equidistantCurve <- function(x,n=NULL,open=TRUE,subsample=0,increment=2,smoothit=0,mesh=NULL,iterations=1) {
    

    
    if (!open)
        x <- rbind(x,x[1,])
    if (is.null(n)) {
        n <- nrow(x)
        if (!open)
            n <- n-1
        }
    up=increment
    k <- nrow(x)
    if (increment > 1 && k < n) {
        while(floor(k*up) <= n) {
            x <- curveSmooth(equidistantCurveHelper(x,floor(k*up),subsample = subsample),iterations=smoothit)
            k <- nrow(x)
        }
        if (k < n)
            x <- curveSmooth(equidistantCurveHelper(x,n),iterations=smoothit)
    } else {
        x <- curveSmooth(equidistantCurveHelper(x,n,subsample = subsample),iterations=smoothit)
    }
    if (iterations > 1) {
        count <- 1
        while(count < iterations) {
            x <- equidistantCurve(x,open=open,subsample=subsample,smoothit=smoothit,mesh=mesh)
            count <- count+1
        }
    }
    
    if (!is.null(mesh))
        x <- vert2points(projRead(x,mesh))
    return(x)
}


equidistantCurveHelper <- function(x,n=100,subsample=0,seed=42) {
    if (subsample > 0 && subsample < nrow(x))
        x <- x[sort(fastKmeans(x,k=subsample)$selected),]
    dists <- x[-1,]-x[-nrow(x),]
    dists <- c(0,sqrt(rowSums(dists^2)))
    meandist <- mean(dists)
    m <- ncol(x)
    k <- nrow(x)
    set.seed(42)
    noise <- rnorm(m*k,sd=meandist/1e3)
    noisen <- rnorm(n*m,sd=meandist/1e3)
    ## print(range(noise))
    flatCurve <- matrix(noise,nrow(x),ncol(x))
    flatCurve[c(1,nrow(x)),] <- 0
    flatCurve[,1] <- cumsum(dists)
    xout <- matrix(noisen,n,ncol(x))
    xout[,1] <- seq(from=0,to=flatCurve[nrow(x),1],length.out = n)
    xout <- tps3d(xout,flatCurve,x,lambda=0)
    return(xout)
}

curveSmooth <- function(x,iterations=1) {
    count <- 0
    xsmooth <- x
    while (count < iterations) {
        k <- nrow(x)
        xcore <- x[-c(1,k),]
        xcoreleft <- x[-c(1:2) ,]
        xcoreright <- x[-c(k,k-1),]
        xsmooth <- (xcore+xcoreleft+xcoreright)/3
        xsmooth <- rbind(x[1,],xsmooth,x[k,])
        count <- count+1
    }
    return(xsmooth)
}


#' sort curvepoints by using the subsequent neighbours
#'
#' sort curvepoints by using the subsequent neighbours
#' @param x k x m matrix containing the 2D or 3D coordinates
#' @param k number of nearest neighbours to look at. Set high for very irregularly clustered curves.
#' @param start integer: which row of x to use as a starting point. If NULL, it is assumed that the curve is open and the point where the angle between the two nearest neighbours is closest will be chosen.
#' @return
#' \item{xsorted}{matrix with coordinates sorted along a curve}
#' \item{index}{vector containing the sorting indices}
#' @examples
#'
#' ## generate a curve from a polynome
#' x <- c(32,64,96,118,126,144,152.5,158)
#' y <- c(99.5,104.8,108.5,100,86,64,35.3,15)
#' fit <- lm(y~poly(x,2,raw=TRUE))
#' xx <- seq(30,160, length=50)
#' layout(matrix(1:3,3,1))
#' curve <- cbind(xx,predict(fit, data.frame(x=xx)))
#' ## permute order
#' set.seed(42)
#' plot(curve);lines(curve)
#' curveunsort <- curve[sample(1:50),]
#' ## now the curve is scrambled
#' plot(curveunsort);lines(curveunsort,col=2)
#' curvesort <- sortCurve(curveunsort)
#' ## after sorting lines are nice again
#' plot(curvesort$xsorted);lines(curvesort$xsorted,col=3)
#' @export
sortCurve <- function(x,k=5, start=NULL) {
    x <- as.matrix(x)
    nn <- vcgKDtree(x,x,k=min(k,nrow(x)))
    if (is.null(start)){
        ## get direction vectors for each point
        dir1 <- x[nn$index[,2],]-x
        dir2 <- x[nn$index[,3],]-x
        angles <- angM(dir1,dir2)
        start <- which.min(angles)
    }
    index0 <- start    
    for (i in 1:(nrow(x)-1)) {
        newind <- nn$index[index0[i],min(which(!nn$index[index0[i],] %in% index0))]
        index0 <- c(index0,newind)
    }
    return(list(xsorted=x[index0,],index=index0))
}



