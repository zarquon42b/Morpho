NNshapeReg <- function(x,y,n,mahalanobis=FALSE)

  {
    out <- y
    for (i in 1:dim(y)[1])
      {
        weighcalc <- proc.weight(x,n,i,mahalanobis=mahalanobis,report=F)$data
        ws <- diag(weighcalc$weight)
        out[i,] <- apply(t(t(y[weighcalc$nr,])%*%ws),2,mean)
      }
    return(out)
  }
    
