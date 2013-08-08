.bindArr2 <- function(x,y,along=1)
    {
        if (is.matrix(x))
            x <- array(x,dim=c(dim(x),1))
        if (is.matrix(y))
            y <- array(y,dim=c(dim(y),1))              
        xdim <- dim(x)
        ydim <- dim(y)
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
        return(newarr)
    }
bindArr <- function(...,along=1)
    {
        args <- list(...)
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
