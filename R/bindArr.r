bindArr <- function(x,y,along=1)
  {
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
              {
               
                newarr[,,i] <- rbind(x[,,i],y[,,i])
              }
            else
              {
                
                newarr[,,i] <- cbind(x[,,i],y[,,i])
              }
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
