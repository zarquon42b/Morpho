cylinder <- function(x,dirs,length)
  {
    circ <- seq(from=0,to=1.9, by=0.1)*pi
    data <- matrix(0,20,3)
    data[,1] <- sin(circ)
    data[,2] <- cos(circ)
   # points3d(data)
   # data2 <- data
   # data2[,3] <- length
    lc <- dim(data)[1]
    it <-NULL
    for (i in 1:(lc-1))
      {
        face0 <- c(i,i+1,i+lc)
        face1 <- c(i+1,i+lc+1,i+lc)
        it <- rbind(it,face0,face1)
      }
    it <- rbind(it,c(lc,1,lc+1))
    it <- rbind(it,c(2*lc,lc,lc+1))
### close lids ###
    for (i in 1:(lc-1))
     {
        it <- rbind(it,c(lc*2+1,i,i+1))
       # face1 <- c(i+1,i+lc+1,i+lc)
       # it <- rbind(it,face0,face1)
      }
    it <- rbind(it,c(lc*2+1,lc,1))
   # it <- rbind(it,face0,face1)
 for (i in (lc+1):(2*lc-1))
     {
        it <- rbind(it,c(lc*2+2,i,i+1))
       # face1 <- c(i+1,i+lc+1,i+lc)
       # it <- rbind(it,face0,face1)
      }
    it <- rbind(it,c(lc*2+2,lc*2,lc+1))

    dirs <- dirs/sqrt(sum(dirs^2))
    yz <- tanplan(dirs)
    rmat <- cbind(yz$z,yz$y,dirs)
    data <- t(rmat%*%t(data))
    data2 <- t(apply(data,1,function(x){x <- x+length*dirs;return(x)}))
    datavb <- rbind(data,data2)
    x0 <- c(0,0,0)
    x1 <- x0+length*dirs
    datavb <- rbind(datavb,x0,x1)
    cylvb <-  t(cbind(datavb,1))
    cyl <- list()
    class(cyl) <- "mesh3d"
    cyl$vb <- cylvb
    cyl$it <- t(it)
    cyl <- conv2backf(cyl)
    cyl <- adnormals(cyl)
    cyl <- translate3d(cyl,x[1],x[2],x[3])
    return(cyl)
  }
