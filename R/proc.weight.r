#' calculate weights inverse to the distances from the specified observation.
#' 
#' for calculation of a shape model by averaging the observations neighbouring
#' the configuration in question, it is necessary to calculate weights by
#' similarity.
#' 
#' distances of zero will get a weight of 1e12 (this is scaled to all weights
#' summing to one), thus weights for observations further away are converging
#' to zero.
#' 
#' @param data array containing landmark configurations
#' @param number integer: how many of the neighbours are to be involved.
#' @param ref integer: position in the array that is used as reference.
#' @param report logical: require report about name of the reference.
#' @param reg numeric: regularise mahalanobis distance by adding reg to the
#' diagonal of eigenvalues of the covariance matrix.
#' @param log logical: use the logarithm of the distances.
#' @param mahalanobis logical: use mahalanobis distance.
#' @return
#' \item{data }{dataframe containing id, procrustes/mahalanobis distance
#' and weight according to the reference}
#' \item{reference }{returns observations' names if available}
#' \item{rho.all }{dataframe containing distances to references of all observations}
#' @examples
#' 
#' library(shapes)
#' proc <- procSym(gorf.dat)
#' ##get weights for the the four specimen closest to the first observation.
#' weights <- proc.weight(proc$rotated,4,1)
#' 
#' ##estimate the first specimen by weighted neighbour shapes.
#' estim <- proc$mshape*0;
#' for (i in 1:4)
#' {estim <-estim+proc$rotated[,,weights$data$nr[i]]*weights$data$weight[i]}
#' 
#' ### visualise
#' plot(estim,asp=1)## show estimation
#' points(proc$rotated[,,1],col=3)##show original
#' 
#' @export
proc.weight <- function(data,number,ref,report=TRUE,reg=0,log=FALSE,mahalanobis=FALSE)
{
  lmdat <- FALSE
  col <- 3
  rho <- 0
  if (length(dim(data))==3)
    {
      data <- vecx(data)
      lmdat <- TRUE
    }
  l <- dim(data)[1]
  obsnames <- dimnames(data)[[1]]
  if (lmdat && !mahalanobis)
    {
      for(i in 1:l)
          rho[i] <- angle.calc(data[ref,],data[i,])
      if (is.null(obsnames))
          id <- as.character(c(1:l))
      else
          id <- obsnames
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
              covtmp <- armaGinv(covtmp)
              cat("singular covariance matrix: using general inverse\n")
            }
        }
      else
          covtmp <- diag(ncol(data))
      rho <- sqrt(mahalanobis(data,cov=covtmp,center=data[ref,],inverted = TRUE))
    }
  
  if (is.null(obsnames))
      id <- as.character(c(1:l))
  else
      id <- obsnames
  if (log)
      rho <- log(rho)
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
      cat(paste("  reference was",id[ref],"\n"))
  
  return(list(data=out,reference=id[ref],rho.all=data))

}
