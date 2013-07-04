NNshapeReg <- function(x,y=NULL, n=3, mahalanobis=FALSE,mc.cores = detectCores())
  {
      if (is.null(y))
          y <- x
      outdim <- dim(y)
      if (length(dim(x)) == 3)
          x <- vecx(x)
      if (length(dim(y)) == 3)
          y <- vecx(y)
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

      if (length(outdim) == 3)
          {
              out1 <- array(NA, dim=outdim)
              for (i in 1:outdim[3])
                  out1[,,i] <- matrix(out[i,],outdim[1],outdim[2])
                   
              out <- out1
          }
    return(out)
  }

