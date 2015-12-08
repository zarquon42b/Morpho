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
#' @param useCor if TRUE, the correlation matrix instead of the covariance matrix is used.
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
#' \item{xlm}{linear model: \code{lm(Xscores ~ Yscores - 1)}}
#' \item{ylm}{linear model: \code{lm(Yscores ~ Xscores - 1)}}
#' @author Stefan Schlager
#' @seealso \code{\link{plsCoVar}, \link{getPLSfromScores}, \link{predictPLSfromScores}, \link{getPLSscores}, \link{predictPLSfromData},\link{svd}}
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
#' layout(1)
#' predPLS <- predictPLSfromData(pls1,y=proc$rotated[5:8,,1])
#' ## show differences between prediction and original
#' deformGrid2d(predPLS,proc$rotated[1:4,,1],pch=19)
#' ##plot the complete first config
#' points(proc$rotated[,,1])
#'
#' ##show effects of first latent variable
#' plsEffects <- plsCoVar(pls1,i=1)
#' deformGrid2d(plsEffects$x[,,1],plsEffects$x[,,2])##show on x
#' deformGrid2d(plsEffects$y[,,1],plsEffects$y[,,2],add=TRUE,pch=19)##show on y
#'
#' ##show effects of 2nd latent variable
#' plsEffects2 <- plsCoVar(pls1,i=2)
#' deformGrid2d(plsEffects2$x[,,1],plsEffects2$x[,,2])##show on x
#' deformGrid2d(plsEffects2$y[,,1],plsEffects2$y[,,2],add=TRUE,pch=19)##show on y
#' @export
pls2B <- function(x, y, tol=1e-12, same.config=FALSE, rounds=0,useCor=FALSE, mc.cores=parallel::detectCores()) {
    
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

        mycorfun <- cov
        if (useCor)
            mycorfun <- cor
        cova <- mycorfun(cbind(x,y))
        
        
        svd.cova <- svd(cova[1:xdim[2],c((xdim[2]+1):(xdim[2]+ydim[2]))])

        svs <- svd.cova$d
        svs <- svs/sum(svs)
        svs <- svs[which(svs > tol)]
        
        covas <- svs*100
        l.covas <- length(covas)
        svd.cova$d <- svd.cova$d[1:l.covas]
        svd.cova$u <- svd.cova$u[,1:l.covas]
        svd.cova$v <- svd.cova$v[,1:l.covas]
        Xscores <- x%*%svd.cova$u#pls scores of x
        Yscores <- y%*%svd.cova$v #pls scores of y
        
        
### calculate correlations between pls scores
        cors <- 0
        for(i in 1:length(covas))
            cors[i] <- cor(Xscores[,i],Yscores[,i])
        

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
                cova.tmp <- mycorfun(cbind(x1[x.sample,],y1[y.sample,]))
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
### get weights
        xlm <- lm(Xscores ~ Yscores -1)
        ylm <- lm(Yscores ~ Xscores -1)
### create covariance table
        Cova <- data.frame(svd.cova$d[1:l.covas],covas,cors,p.values)
        colnames(Cova) <- c("singular value","% total covar.","Corr. coefficient", "p-value")
        out <- list(svd=svd.cova,Xscores=Xscores,Yscores=Yscores,CoVar=Cova)
        out$x <- xorig
        out$y <- yorig
        out$xcenter <- attributes(x)$"scaled:center"
        out$ycenter <- attributes(y)$"scaled:center"
        out$xlm <- xlm
        out$ylm <- ylm
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

#' compute changes associated with 2-Block PLS-scores 
#'
#' compute changes associated with 2-Block PLS-scores 
#'
#' @param pls output of pls2B
#' @param x scores associated with dataset x in original pls2B
#' @param y scores associated with dataset y in original pls2B
#' @return returns data in the original space associated with the specified values.
#' @details other than \code{\link{predictPLSfromScores}}, providing Xscores will not compute predictions of y, but the changes in the original data \code{x} that is associated with the specific scores
#' @export
getPLSfromScores <- function(pls,x,y) {
    if (!missing(x) && !missing(y))
        stop("either x or y must be missing")
        svdpls <- pls$svd

    if (missing(y)) {
        if (is.vector(x) || length(x) == 1) {
            xl <- length(x)
            x <- t(x)
        }
        else if (is.matrix(x))
            xl <- ncol(x)

        out <- t(svdpls$u[,1:xl]%*%t(x))
        out <- sweep(out,2,-pls$xcenter)
        if (length(dim(pls$x)) == 3) {
            if (is.matrix(x) && nrow(x) > 1) {
                out <- vecx(out,revert = T,lmdim = dim(pls$x)[2])
                #dimnames(out) <- dimnames(pls$y)
            } else {
                out <- matrix(out,dim(pls$x)[1],dim(pls$x)[2])
            }
        }
        return(out)
    }
    if (missing(x)) {
        if (is.vector(y) || length(y) == 1) {
            xl <- length(y)
            y <- t(y)
        } else if (is.matrix(y))
            xl <- ncol(y)
        
        out <- t(svdpls$v[,1:xl]%*%t(y))
        out <- sweep(out,2,-pls$ycenter)
        if (length(dim(pls$y)) == 3) {
            if (is.matrix(y) && nrow(y) > 1) {
                out <- vecx(out,revert = T,lmdim = dim(pls$y)[2])
                #dimnames(out) <- dimnames(pls$y)
            } else {
                out <- matrix(out,dim(pls$y)[1],dim(pls$y)[2])
            }
        }
        return(out)
    }
    
}
    

    

#' predict data from 2-Block PLS-scores
#'
#' predict data from 2-Block PLS-scores
#' @param pls output of pls2B
#' @param x scores associated with dataset x in original pls2B
#' @param y scores associated with dataset y in original pls2B
#' @note either x or y must be missing. If x-scores are provided, the yscores will be estimated and the predictions calculated.
#' @return returns an array/matrix of landmarks or original values, depending on input for computing \code{pls}
#' @seealso \code{\link{pls2B}, \link{getPLSscores},\link{predictPLSfromData}}
#' 
#' @export
predictPLSfromScores <- function(pls,x,y) {
    if (!missing(x) && !missing(y))
        stop("either x or y must be missing")
    
    svdpls <- pls$svd
    if (missing(y)) {
        if (is.vector(x) || length(x) == 1) {
            xl <- length(x)
            x <- t(x)
        }
        else if (is.matrix(x))
            xl <- ncol(x)

        if (xl == 1)
            yest <- pls$ylm$coefficients[1:xl,]*x
        else
            yest <- t(t(pls$ylm$coefficients[1:xl,])%*%t(x))
        out <- t(svdpls$v%*%t(yest))
        out <- sweep(out,2,-pls$ycenter)
        if (length(dim(pls$y)) == 3) {
            if (is.matrix(x) && nrow(x) > 1) {
                out <- vecx(out,revert = T,lmdim = dim(pls$x)[2])
                #dimnames(out) <- dimnames(pls$y)
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
        if (xl == 1)
            xest <- pls$xlm$coefficients[1:xl,]*y
        else
            xest <- t(t(pls$xlm$coefficients[c(1:xl),])%*%t(y))
        out <- t(svdpls$u%*%t(xest))
        
        out <- sweep(out,2,-pls$xcenter)
        if (length(dim(pls$x)) == 3) {
            if (is.matrix(y) && nrow(y) > 1) {
                out <- vecx(out,revert = T,lmdim = dim(pls$x)[2])
                #dimnames(out) <- dimnames(pls$x)
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
#' @param ncomp number of (latent) components to use for prediction.
#' @note either x or y must be missing
#' @return returns an array/matrix/vector of predictions - depending on input for computing \code{pls}
#' @seealso \code{\link{pls2B}, \link{getPLSscores},\link{predictPLSfromScores}}
#' @examples
#' ##see examples in pls2B
#' @export
predictPLSfromData <- function(pls,x,y,ncomp=NULL) {
    if (!missing(x) && !missing(y))
        stop("either x or y must be missing")
    if (is.null(ncomp))
        ncomp <- ncol(pls$Xscores)

    if (missing(y)) {
        scores <- getPLSscores(pls,x=x)[,1:ncomp]
        out <- predictPLSfromScores(pls,x=scores)

    }
    if (missing(x)) {
        scores <- getPLSscores(pls,y=y)[,1:ncomp]
        out <- predictPLSfromScores(pls,y=scores)
    }
    return(out)
}


#' Get the shape changes from pls2B associated with each latent variable
#'
#' Get the shape changes from pls2B associated with each latent variable
#' @param pls output of pls2B
#' @param i integer: which latent variable to show. E.g. i=3 will show the changes associated with the 3rd latent variable.
#' @param sdx standard deviation on the xscores. sdx=3 will show the effecs of +3sd vs -3sd
#' @param sdy standard deviation on the yscores. sdy=3 will show the effecs of +3sd vs -3sd
#' @return \item{x}{matrix/array with reconstructed x}
#' @return \item{y}{matrix/array with reconstructed y, with each prediction named accordingly: e.g. neg_x_sd_3 means the prediction of x at a score of \code{-3*sd(Xscores)}}. 
#' @seealso \code{\link{pls2B}, \link{getPLSfromScores}, \link{predictPLSfromScores}, \link{getPLSscores}, \link{predictPLSfromData},\link{svd}}
#' @export 
plsCoVar <- function(pls,i,sdx=3,sdy=3) {
    
    x <- t(t(c(-1,1)*sdx*sd(pls$Xscores[,i])))
    y <- t(t(c(-1,1)*sdy*sd(pls$Yscores[,i])))
    x0 <- matrix(0,2,i); x0[,i] <- x
    y0 <- matrix(0,2,i); y0[,i] <- y
    xnames <-  paste(c("neg","pos"),"x_sd",sdx,sep="_")
    ynames <-  paste(c("neg","pos"),"y_sd",sdy,sep="_")
    pls1x <- getPLSfromScores(pls,x=x0)
    if (is.matrix(pls1x))
        rownames(pls1x) <- xnames
    else
        dimnames(pls1x)[[3]] <- xnames
    pls1y <- getPLSfromScores(pls,y=y0)
    if (is.matrix(pls1y))
        rownames(pls1y) <- ynames
    else
        dimnames(pls1y)[[3]] <- ynames

    pls1out <- list(x=pls1x,y=pls1y)
    return(pls1out)
}
