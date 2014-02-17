#' correlation between shape space and tangent space
#' 
#' performs a partial Procrustes superimposition of landmark data and
#' calculates the correlation between tangent and shape space.
#' 
#' 
#' @param dataarray Input k x m x n real array, where k is the number of
#' points, m is the number of dimensions, and n is the sample size.
#' @param plot Logical: whether to plot the distances between observations.
#' @param main character string: Title of the plot.
#' @param rho chose how to calculate distances in shape space. Options:
#' "riemdist"=Riemannian distance (function from the shapes package-takes along
#' time to calculate), "angle"=calculates the angle between shape vectors,
#' "sindist"=sinus of length of residual vector between shape vectors.
#' @param dist.mat.out Logical: If TRUE, output will contain distance matrices.
#' @return
#' \item{cor }{correlation coefficient between distances in shape space and
#' tangent space}
#' \item{procSS }{Procrustes Sums of Squares (of full procrustes distance)}
#' \item{tanSS }{Tangent Sums of Squares}
#' \item{rhoSS }{Procrustes Sums of Squares (of angle)}
#' \item{euc.dist }{distance matrix of euclidean distance in Tangent space}
#' \item{proc.dist }{distance matrix of Procrustes distance in Shape space}
#' @author Stefan Schlager
#' @seealso \code{\link{regdist}}
#' 
#' @examples
#' 
#' library(shapes)
#' regdist(gorf.dat)
#' 
#' @export
regdist <- regdist.raw <- function(dataarray, plot=TRUE, main="", rho="angle", dist.mat.out=FALSE)
{     proc <- procSym(dataarray,scale=FALSE)
      x <- proc$rotated
      n <- dim(x)[3]
      m <- dim(x)[2]
      k <- dim(x)[1]
      y <- proc$orpdata

      qm <- dist(t(matrix(x,k*m,n)))  #calc  dist. between rotated config
      procdis <- sum(qm^2)/n
      procdistmat <- matrix(NA,n,n) #calc rho from angle between rotated configs
      for (i in 1:n)
          for (j in 1:n) {
              if (rho=="riemdist") {
                  procdistmat[i,j] <- kendalldist(x[,,i],x[,,j])  # riemann dist.
              } else if (rho=="angle") 
                  procdistmat[i,j] <- angle.calc(x[,,i],x[,,j])
          }
      if (rho == "sindist")
          procvec <- asin(qm)
      else
          procvec <- as.dist(procdistmat)
      
      procdis2 <- sum(procvec^2)/n
      em <- dist(t(matrix(y,k*m,n)))
      euvec <- (em)
      eudis <- sum(euvec^2)/n
      correlation <- cor(euvec,procvec)^2
      
      if (plot==TRUE)
          plot(euvec,procvec,asp=1,xlab="euclid. dist. in tangentspace",ylab=paste("rho as",rho),main=main)
      abline(0,1,col="grey50")
      
      if (dist.mat.out)
          return(list(cor=correlation,procSS=procdis,tanSS=eudis,rhoSS=procdis2,euc.dist=em,proc.dist=procvec))
      else
          return(list(cor=correlation,procSS=procdis,tanSS=eudis,rhoSS=procdis2))
  }
