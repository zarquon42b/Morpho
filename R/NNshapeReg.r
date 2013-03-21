NNshapeReg <- function(x,y,n,mahalanobis=FALSE,mc.cores = detectCores())
  {
    i <- NULL
    win <- FALSE
     if(.Platform$OS.type == "windows")
      win <- TRUE
    else
      registerDoParallel(cores=mc.cores)### register parallel backend
    out <- y
    estfun <- function(i)
      {
        weighcalc <- proc.weight(x,n,i,mahalanobis=mahalanobis,report=F)$data
        ws <- diag(weighcalc$weight)
        tmpres <- apply(t(t(y[weighcalc$nr,])%*%ws),2,sum)
        return(tmpres)
      }
    if (win)
      out <- foreach(i=1:dim(x)[1],.combine=rbind) %do% estfun(i)
    else
      out <- foreach(i=1:dim(x)[1],.combine=rbind) %dopar% estfun(i)
    return(out)
  }

