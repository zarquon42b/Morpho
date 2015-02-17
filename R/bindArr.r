.bindArr2 <- function(x,y,along=1)
    {
        if (is.matrix(x))
            x <- array(x,dim=c(dim(x),1))
        if (is.matrix(y))
            y <- array(y,dim=c(dim(y),1))

        if (length(y) == 1 && along %in% c(1,2)) {
            if (along == 2)
                y <- array(y,dim=c(dim(x)[2],1,dim(x)[3]))
            else
                y <- array(y,dim=c(1,dim(x)[2],dim(x)[3]))
        }
                                 
        xdim <- dim(x)
        ydim <- dim(y)
        outnames <- list()
        xnames <- dimnames(x)
        ynames <- dimnames(y)
        for (i in (1:3)) {
            check <- is.null(xnames[[i]])
            check[2] <- is.null(ynames[[i]])
            if (!prod(check)) {
                if (i != along) {
                    tmpsep <- "_"
                    if (sum(check)) {
                        tmpsep=""
                    } else {
                        if (prod(xnames[[i]] == ynames[[i]])) {
                            xnames[[i]] <- ""
                            tmpsep=""
                        }
                    }
                    outnames[[i]] <- paste(xnames[[i]],ynames[[i]],sep=tmpsep)
                } else {
                    if (!prod(check)) {
                        if (check[1])
                            outnames[[along]] <- c(paste0("X",1:xdim[along]),ynames[[along]])
                        else if (check[2])
                            outnames[[along]] <- c(xnames[[along]],paste0("X",1:ydim[along]))
                        else
                            outnames[[along]] <- append(xnames[[along]],ynames[[along]])
                    }
                }
            } else {
                outnames[[i]] <- NULL
            }
        }            
        newalong <- xdim[along]+ydim[along]
        if (along %in% 1:2)
            {
                if (along == 1)
                    newarr <- array(NA,c(newalong,xdim[2:3]))
                else
                    newarr <- array(NA,c(xdim[1],newalong,xdim[3]))
                for(i in 1:xdim[3])
                    {
                        if (along==1)
                            newarr[,,i] <- rbind(x[,,i],y[,,i])
                        else
                            newarr[,,i] <- cbind(x[,,i],y[,,i])
                    }
            }
        else
            {
                newarr <- array(NA,c(xdim[1:2],newalong))
                newarr[,,1:xdim[3]] <- x
                newarr[,,(xdim[3]+1):newalong] <- y
            }
        
        dimnames(newarr) <- outnames
        return(newarr)
    }


#' concatenate multiple arrays/matrices
#' 
#' concatenate multiple 3-dimensional arrays and/or 2-dimensional matrices to
#' one big array
#' 
#' 
#' @param \dots matrices and/or arrays with appropriate dimensionality to
#' combine to one array, or a single list containing suitable matrices, or arrays).
#' @param along dimension along which to concatenate.
#' @details dimnames, if present and if differing between entries, will be concatenated, separated by a "_".
#' @return returns array of combined matrices/arrays
#' @seealso \code{\link{cbind}}, \code{\link{rbind}}, \code{\link{array}}
#' 
#' @examples
#' 
#' A <- matrix(rnorm(18),6,3)
#' B <- matrix(rnorm(18),6,3)
#' C <- matrix(rnorm(18),6,3)
#' 
#' #combine to 3D-array
#' newArr <- bindArr(A,B,C,along=3)
#' #combine along first dimension
#' newArr2 <- bindArr(newArr,newArr,along=1)
#' 
#' 
#' 
#' @export
bindArr <- function(...,along=1)
    {
        args <- list(...)
        if (length(args) == 1 && is.list(args[[1]]))
             args <- (...)
       
        argc <- length(args)
        if (argc < 2)
            stop("at least two arguments needed")
        newarr <- .bindArr2(args[[1]],args[[2]], along=along)
        if (argc > 2) {
            for (i in 3:argc)
                newarr <- .bindArr2(newarr, args[[i]],along=along)
        }
        return(newarr)
    }
