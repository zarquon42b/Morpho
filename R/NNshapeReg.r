NNshapeReg <- function(x,y,n,mahalanobis=FALSE,mc.cores = detectCores())

  {
    i <- NULL
    registerDoParallel(mc.cores)
    out <- y
    estfun <- function(i)
      {
        weighcalc <- proc.weight(x,n,i,mahalanobis=mahalanobis,report=F)$data
        ws <- diag(weighcalc$weight)
        tmpres <- apply(t(t(y[weighcalc$nr,])%*%ws),2,sum)
        return(tmpres)
      }
    out <- foreach(i=1:dim(x)[1],.combine=rbind) %dopar% estfun(i)
    return(out)
  }
    
