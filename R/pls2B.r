#' Two-Block partial least square regression.
#' 
#' Performs a Two-Block PLS on two sets of data and assesses the significance
#' of each score by permutation testing
#' 
#' The Two-Block PLS tries to find those linear combinations in each block
#' maximising the covariance between blocks. The significance of each linear
#' combination is assessed by comparing the singular value to those obtained
#' from permuted blocks. If both blocks contain landmarks superimposed
#' TOGETHER, the option \code{same.config=TRUE} requests superimposition of the
#' permuted configurations (i.e. where the the landmarks of block \code{x} are
#' replaced by corresponding landmarks of other specimen.
#' 
#' @param y array containing superimposed landmark data of the first block.
#' Matrices are also allowed but the option 'same.config' will not work.
#' @param x array containing superimposed landmark data second block.Matrices
#' are also allowed but the option 'same.config' will not work.
#' @param tol threshold for discarding singular values.
#' @param same.config logical: if \code{TRUE} each permutation includes new
#' superimposition of permuted landmarks. This is necessary if both blocks
#' originate from landmarks that are superimposed together.
#' @param rounds rounds of permutation testing.
#' @param mc.cores integer: determines how many cores to use for the
#' computation. The default is autodetect. But in case, it doesn't work as
#' expected cores can be set manually. Parallel processing is disabled on
#' Windows due to occasional errors.
#' @return
#' \item{svd }{singular value decomposition (see \code{\link{svd}}) of the
#' 'common' covariance block}
#' \item{Xscores }{PLS-scores of x}
#' \item{Yscores }{PLS-scores of y}
#' \item{CoVar }{Dataframe containing singular values, explained
#' covariation, correlation coeffictient between PLS-scores and p-values}
#' @author Stefan Schlager
#' @seealso \code{\link{svd}}
#' @references Rohlf FJ, Corti M. 2000. Use of two-block partial least-squares
#' to study covariation in shape. Systematic Biology 49:740-753.
#' @examples
#' 
#' library(shapes)
#' ### very arbitrary test:
#' ### check if first 4 landmarks covaries with the second 4
#' proc <- procSym(gorf.dat)
#' ## we do only 50 rounds to minimize computation time
#' \dontrun{#same.config takes too long for CRAN check
#' pls1 <- pls2B(proc$rotated[1:4,,],proc$rotated[5:8,,],
#'               same.config=TRUE,rounds=50,mc.cores=2)
#' }
#' pls1 <- pls2B(proc$rotated[1:4,,],proc$rotated[5:8,,],
#'               same.config=FALSE,rounds=50,mc.cores=1)
#' pls1
#' layout(matrix(1:4,2,2,byrow=TRUE))
#' for(i in 1:4)
#'  plot(pls1$Xscores[,i]~pls1$Yscores[,i])
#'
#'
#' ## predict first 4 landmarks from second 4 for first config
#' 
#' predPLS <- predictPLSfromData(pls1,y=proc$rotated[5:8,,1])
#' ## show differences between prediction and original
#' deformGrid2d(predPLS,proc$rotated[1:4,,1],pch=19)
#' ##plot the complete first config
#' points(proc$rotated[,,1])
#' 
#' @export
pls2B <- function(x, y, tol=1e-12, same.config=FALSE, rounds=0, mc.cores=parallel::detectCores(),scale=FALSE) {
    
        landmarks <- FALSE
        xorig <- x
        yorig <- y
        win <- FALSE
        if(.Platform$OS.type == "windows")
            win <- TRUE
        else
            registerDoParallel(cores=mc.cores)### register parallel backend
        
        if (length(dim(x)) == 3) {
            landmarks <- TRUE
            x <- vecx(x)
        }
        if (length(dim(y)) == 3)
            y <- vecx(y)
        else
            landmarks <- FALSE
        
        xdim <- dim(x)
        ydim <- dim(y)

        if (same.config && !landmarks)
            warning("the option same.config requires landmark array as input")
        
        x <- scale(x,scale = F)
        y <- scale(y,scale = F)
        if (!scale)
            cova <- cov(cbind(x,y))
        else
            cova <- cor(cbind(x,y))
        
        svd.cova <- svd(cova[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])

        svs <- svd.cova$d
        svs <- svs/sum(svs)
        svs <- svs[which(svs > tol)]
        
        covas <- svs*100
        l.covas <- length(covas)
        svd.cova$d <- svd.cova$d[1:l.covas]
        svd.cova$u <- svd.cova$u[,1:l.covas]
        svd.cova$v <- svd.cova$v[,1:l.covas]
        z1 <- x%*%svd.cova$u#pls scores of x
        z2 <- y%*%svd.cova$v #pls scores of y
        
        
### calculate correlations between pls scores
        cors <- 0
        for(i in 1:length(covas))
            cors[i] <- cor(z1[,i],z2[,i])
        

### Permutation testing
        permupls <- function(i)
            {
                x.sample <- sample(1:xdim[1])
                y.sample <- sample(x.sample)
                if (same.config && landmarks) {
                    tmparr <- .bindArr2(xorig[,,x.sample],yorig[,,y.sample],along=1)
                    tmpproc <- ProcGPA(tmparr,silent=TRUE)
                    x1 <- vecx(tmpproc$rotated[1:dim(xorig)[1],,])
                    y1 <- vecx(tmpproc$rotated[1:dim(yorig)[1],,])
                } else {
                    x1 <- x
                    y1 <- y
                }
                x1 <- scale(x1,scale = F)
                y1 <- scale(y1,scale = F)
                cova.tmp <- cov(cbind(x1[x.sample,],y1[y.sample,]))
                svd.cova.tmp <- svd(cova.tmp[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])
                svs.tmp <- svd.cova.tmp$d
                return(svs.tmp[1:l.covas])
            }
        p.values <- rep(NA,l.covas)
        if (rounds > 0) {
            if (win)
                permuscores <- foreach(i = 1:rounds, .combine = cbind) %do% permupls(i)
            else
                permuscores <- foreach(i = 1:rounds, .combine = cbind) %dopar% permupls(i)
            
            p.val <- function(x,rand.x)
                {
                    p.value <- length(which(rand.x >= x))
                    
                    if (p.value > 0)
                        p.value <- p.value/rounds
                    else
                        p.value <- 1/rounds
                    gc()
                    return(p.value)
                }
            
            for (i in 1:l.covas)
                p.values[i] <- p.val(svd.cova$d[i],permuscores[i,])
            
        }
### create covariance table
        Cova <- data.frame(svd.cova$d[1:l.covas],covas,cors,p.values)
        colnames(Cova) <- c("singular value","% total covar.","Corr. coefficient", "p-value")
        out <- list(svd=svd.cova,Xscores=z1,Yscores=z2,CoVar=Cova)
        out$x <- xorig
        out$y <- yorig
        out$xcenter <- attributes(x)$"scaled:center"
        out$ycenter <- attributes(y)$"scaled:center"
        class(out) <- "pls2B"
        return(out)
    }

#' @export
print.pls2B <- function(x,...) {
    cat("  Covariance explained by the singular values\n\n")
    df <- x$CoVar
    df <- df[,colSums(is.na(df)) != nrow(df)]
    print( df,row.names=FALSE)
}

#' predict data from 2-Block PLS-scores
#'
#' predict data from 2-Block PLS-scores
#' @param pls output of pls2B
#' @param x scores associated with dataset x in original pls2B
#' @param y scores associated with dataset x in original pls2B
#' @note either x or y must be missing
#' @return returns an array/matrix of landmarks or original values, depending on input for computing \code{pls}
#' @seealso \code{\link{pls2B}, \link{getPLSscores},\link{predictPLSfromData}}
#' 
#' @export
predictPLSfromScores <- function(pls,x,y) {
    if (!missing(x) && !missing(y))
        stop("either x or y must be missing")
    
    scalevec <- apply(pls$Yscores,2,sd)/apply(pls$Xscores,2,sd)
    svdpls <- pls$svd
    if (missing(y)) {
        if (is.vector(x) || length(x) == 1) {
            xl <- length(x)
            x <- t(x)
        }
        else if (is.matrix(x))
            xl <- ncol(x)
        lmpred <- lm(pls$Yscores ~ pls$Xscores)
        yest <- predict(lmpred,newdata = list("pls$Xscores"=x))
        #scaledv <- t(scalevec[1:xl]*t(svdpls$v[,1:xl]))
        scaledv <- svdpls$v[,1:xl]
        out <- t(scaledv%*%t(yest))
        out <- sweep(out,2,-pls$ycenter)
        if (length(dim(pls$y)) == 3) {
            if (is.matrix(x) && nrow(x) > 1) {
                out <- vecx(out,revert = T,lmdim = dim(pls$x)[2])
                dimnames(out) <- dimnames(pls$y)
            } else {
                out <- matrix(out,dim(pls$y)[1],dim(pls$y)[2])
            }
        }
        
    }
    
    if (missing(x)) {
        if (is.vector(y) || length(y) == 1) {
            xl <- length(y)
            y <- t(y)
        }
        else if (is.matrix(y))
            xl <- ncol(y)
        lmpred <- lm(pls$Xscores ~ pls$Yscores)
        xest <- predict(lmpred,newdata = list("pls$Yscores"=y))
        #scaledv <- t(scalevec[1:xl]*t(svdpls$v[,1:xl]))
        scaledu <- svdpls$u[,1:xl]
        out <- t(scaledu%*%t(xest))
        out <- sweep(out,2,-pls$xcenter)
        if (length(dim(pls$x)) == 3) {
            if (is.matrix(y) && nrow(y) > 1) {
                out <- vecx(out,revert = T,lmdim = dim(pls$x)[2])
                dimnames(out) <- dimnames(pls$y)
            } else {
                out <- matrix(out,dim(pls$x)[1],dim(pls$x)[2])
            }
        }
        
    }
    return(out)
}


#' compute 2-Block PLS scores for new data 
#'
#' compute 2-Block PLS scores for new data from an existing pls2B
#' @param pls output of pls2B
#' @param x matrix or vector representing new dataset(s) -  same kind as in original pls2B
#' @param y matrix or vector representing new dataset(s) - same kind as in original pls2B
#' @note either x or y must be missing
#' 
#' @return returns a vector of pls-scores
#' @seealso \code{\link{pls2B}, \link{predictPLSfromScores},\link{predictPLSfromData}}
#' @export
getPLSscores <- function(pls,x,y) {
    if (!missing(x) && !missing(y))
        stop("either x or y must be missing")
    
    
    if (missing(y)) {
        if (length(dim(x)) == 3 || (is.matrix(x) && is.matrix(pls$x))) {
            if (length(dim(x)) == 3)
                x <- vecx(x)
            out <- NULL
            for(i in 1:nrow(x))
                out <- rbind(out,getPLSscores(pls,x=x[i,]))
            
        } else {
            if (is.matrix(x))
                x <- as.vector(x)
            x <- x-pls$xcenter
            out <- t(x)%*%pls$svd$u
        }
        
    }
    if (missing(x)) {
        if (length(dim(y)) == 3 || (is.matrix(y) && is.matrix(pls$y))) {
            if (length(dim(y)) == 3)
                y <- vecx(y)
            out <- NULL
            for(i in 1:nrow(y))
                out <- rbind(out,getPLSscores(pls,y=y[i,]))
            
        } else {
            if (is.matrix(y))
                y <- as.vector(y)
            y <- y-pls$ycenter
            out <- t(y)%*%pls$svd$v
            
        }
    }
    return(out)
}
#' predict 2 Block-PLS from new data
#'
#' predict 2 Block-PLS from new data
#' @param pls output of pls2B
#' @param x data in the same format as in original pls2B (for landmarks this can be an array or a matrix and for other data a matrix of a vector)
#' @param y data in the same format as in original pls2B (for landmarks this can be an array or a matrix and for other data a matrix of a vector)
#' @note either x or y must be missing
#' @return returns an array/matrix/vector of predictions - depending on input for computing \code{pls}
#' @seealso \code{\link{pls2B}, \link{getPLSscores},\link{predictPLSfromScores}}
#' @examples
#' ##see examples in pls2B
#' @export
predictPLSfromData <- function(pls,x,y) {
    if (!missing(x) && !missing(y))
        stop("either x or y must be missing")
    

    if (missing(y)) {
        scores <- getPLSscores(pls,x=x)
        out <- predictPLSfromScores(pls,x=scores)

    }
    if (missing(x)) {
        scores <- getPLSscores(pls,y=y)
        out <- predictPLSfromScores(pls,y=scores)
    }
    return(out)
}
