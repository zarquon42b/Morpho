proc.weight<-function(data,number,ref,report=TRUE,reg=0,log=FALSE,mahalanobis=FALSE)
{
  lmdat <- FALSE
  col <- 3
  rho <- 0
  if (length(dim(data))==3)
    {
      data <- vecx(data)
      lmdat <- TRUE
    }
  l<-dim(data)[1]
  obsnames <- dimnames(data)[[1]]
  if (lmdat && !mahalanobis)
    {
      for(i in 1:l){rho[i]<-angle.calc(data[ref,],data[i,])$rho}
      if (is.null(obsnames))
        {id <- as.character(c(1:l))
       }
      else
        {id <- obsnames}
    }
  else
    {
      if (mahalanobis)
        {
          covtmp <- cov(data)
          if (reg != 0)
            {
              {
                eig <- eigen(covtmp,symmetric=TRUE)
                covtmp <- t(eig$vectors)%*%diag(eig$values+reg)%*%eig$vectors
              }
            }
          checksing <- try(covtmp <- solve(covtmp),silent = TRUE)
          if (class(checksing)=="try-error")
            {
              covtmp <- ginv(covtmp)
              cat("singular covariance matrix: using general inverse\n")
            }
        }
      else
        {
          covtmp <- diag(ncol(data))
        }
      rho <- sqrt(mahalanobis(data,cov=covtmp,center=data[ref,],inverted = TRUE))
    }
  
  if (is.null(obsnames))
    {id <- as.character(c(1:l))
   }
  else
    {id <- obsnames}
  if (log)
    {rho <- log(rho)
   }
  nr <- c(1:l)
  data <- data.frame(nr,id,rho)
  dat.sort.i <- data[order(data[,col]),]
  dat.which <- dat.sort.i[2:(number+1),]
  weight <- dat.which[,col]
  if (0 %in% weight)
    weight[which(weight == 0)] <- 1/1e12
  
  weight <- 1/weight
  weight <- weight/sum(weight)
  out <- data.frame(dat.which,weight)
  if (report)
    {
      cat(paste("  reference was",id[ref],"\n"))
    }
  return(list(data=out,reference=id[ref],rho.all=data))

}
