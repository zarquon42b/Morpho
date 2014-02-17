#' estimate missing landmarks
#' 
#' Missing landmarks are estimated by deforming a sample average or a weighted
#' estimate of the configurations most similar onto the deficient
#' configuration. The deformation is performed by a Thin-plate-spline
#' interpolation calculated by the available landmarks.
#' 
#' This function tries to estimate missing landmark data by mapping weighted
#' averages from complete datasets onto the missing specimen. The weights are
#' the inverted Procrustes (see \code{\link{proc.weight}}) distances between
#' the 'comp' closest specimen (using the available landmark configuration).
#' 
#' @param data array containing landmark data
#' @param comp integer: select how many of the closest observations are to be
#' taken to calculate an initial estimate.
#' @param weight logical: requests the calculation of an estimate based on the
#' procrustes distance. Otherwise the sample's consensus is used as reference.
#' @return
#' \item{out }{array containing all data, including fixed configurations - same order as input}
#' \item{mshape }{meanshape - calculated from complete datasets}
#' \item{checklist }{list containing information about missing landmarks}
#' \item{check }{vector containing position of observations in data where at least one missing coordinate was found}
#' @note Be aware that these estimates might be grossly wrong when the missing
#' landmark is quite far off the rest of the landmarks (due to the radial basis
#' function used in the Thin-plate spline interpolation.
#' @author Stefan Schlager
#' @seealso \code{\link{proc.weight}}, \code{\link{tps3d}}
#' @references Bookstein FL. 1989. Principal Warps: Thin-plate splines and the
#' decomposition of deformations IEEE Transactions on pattern analysis and
#' machine intelligence 11.
#' 
#' @examples
#' 
#' require(rgl)
#' require(shapes)
#' data <- gorf.dat
#' ### set first landmark of first specimen to NA
#' data[1,,1] <- NA
#' repair <- fixLMtps(data,comp=5)
#' ### view difference between estimated and actual landmark
#' plot(repair$out[,,1],asp=1,pch=21,cex=0.7,col=2)#estimated landmark
#' points(gorf.dat[,,1],col=3,pch=20)#actual landmark
#' 
#' ## 3D-example:
#' data(boneData)
#' data <- boneLM
#' ### set first and 5th landmark of first specimen to NA
#' data[c(1,5),,1] <- NA
#' repair <- fixLMtps(data,comp=10)
#' ##  view difference between estimated and actual landmark
#' \dontrun{
#' deformGrid3d(repair$out[,,1], boneLM[,,1],ngrid=0)
#' }
#' @export
fixLMtps <- function(data,comp=3,weight=TRUE)
{
  n <- dim(data)[3]
  k <- dim(data)[1]
  m <- dim(data)[2]
  checklist <- list()
  checkvec <- rep(0,n)
  out <- data
  ## check for missing landmarks ###
  for (i in 1:n) {
      count <- 0
      checklist[[i]] <- NA
      for (j in 1:k) {
          if (sum(is.na(data[j,,i]))) {
              count <- count+1
              checklist[[i]][count] <- j
              checkvec[i] <- 1
          }
      }
  }
  ## calc mean of complete configs ###
  
  check <- which(checkvec==1)
  if (!length(check)) {
    message("nothing to fix")
    return(data)
  }
    
  data.c <- data[,,-check]
  ## check if there are enough configs to use weighting option
  if (length(dim(data.c)) < 3) {
      if (is.matrix(data.c)) {
          data.c <- array(data.c,dim=c(dim(data.c),1))
          ngood <- 1
      } else
          stop("there is no complete configuration to use")
  } else {
      ngood <- dim(data.c)[3]
  }
  if (ngood < comp) {
      if (ngood == 0) {
          stop("no complete configuration found")
      } else {
          comp <- ngood
          if (weight)
              warning(paste("only",ngood,"configuration(s) found. comp is set to",ngood,"\n"))
      }
  }
  
  ## rotate incomplete data onto mean ###
  if (ngood > 1) {
      proc.c <- ProcGPA(data.c,silent = TRUE)
      mean0 <- proc.c$mshape
  } else
      mean0 <- data.c[,,1]
  for (i in 1:length(check)) {
      miss <- checklist[[check[i]]]
      if (weight && ngood > 1) { ### calculate weighted estimates of missing data ###
          ## rotate incomplete data onto mean ###
          rotmiss <- rotonto(mean0[-miss,],data[-miss,,check[i]],scale=TRUE)$yrot
          allrot <- bindArr(rotmiss,proc.c$rotated[-miss,,], along=3)
          ## calculate weights according to procrustes distance ###			
          wcalc <- proc.weight(allrot,comp,1,report=FALSE)
          lms <- proc.c$rotated[,,wcalc$data$nr-1]
          if (is.matrix(lms))
              lms <- array(lms,dim=c(dim(lms),1))
          lm.est <- matrix(0,dim(data)[1],m)
          
          for (j in 1:comp) {
              lm.est <- lm.est+lms[,,j]*wcalc$data$weight[j]
          }
          tpsout <- tps3d(lm.est,lm.est[-miss,],data[-miss,,check[i]])
      } else {
          tpsout <- tps3d(mean0,mean0[-miss,],data[-miss,,check[i]])
      }
      out[,,check[i]] <- tpsout
  }
  return(list(out=out,mshape=mean0,checklist=checklist,check=check))
}
