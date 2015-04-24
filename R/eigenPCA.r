eigenPCA <- function(data,tol=.Machine$double.eps)
  {
    center1 <- colMeans(data,2,mean)
    datascale <- scale(data,2,scale=FALSE)
    eigData <- eigen(cov(datascale),symmetric=TRUE)
    eigData$values <- eigData$values[which(eigData$values > tol)] 
    sdev <- sqrt(eigData$values)
    rotation <- eigData$vectors[,1:length(sdev)]
    x <- t(t(rotation)%*%t(datascale))
    return(list(x=x,rotation=rotation,sdev=sdev,center=center1))
  }
