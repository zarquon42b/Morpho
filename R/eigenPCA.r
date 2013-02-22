eigenPCA <- function(data)
  {
    center1 <- apply(data,2,mean)
    datascale <- apply(data,2,scale,scale=FALSE)
    eigData <- eigen(cov(datascale),symmetric=TRUE)
    rotation <- eigData$vectors
    x <- t(t(rotation)%*%t(datascale))
    eigData$values[which(eigData$values < 0)] <- 0
    sdev <- sqrt(eigData$values)
    return(list(x=x,rotation=rotation,sdev=sdev,center=center1))
  }
