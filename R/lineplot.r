#' plot lines between landmarks
#' 
#' add lines connecting landmarks to visualise a sort of wireframe
#' 
#' 
#' @param x matrix containing 2D or 3D landmarks 
#' @param point vector or list of vectors containing rowindices of x,
#' determining which landmarks to connect.
#' @param col color of lines
#' @param lwd line width
#' @param lty line type (only for 2D case)
#' @param line_antialias logical: smooth lines
#' @param add logical: add to existing plot
#' @note works with 2D and 3D configurations
#' @author Stefan Schlager
#' @seealso \code{\link{pcaplot3d}}
#' 
#' @examples
#' 
#' 
#' if (require(shapes)) {
#' ##2D example
#' plot(gorf.dat[,,1],asp=1)
#' lineplot(gorf.dat[,,1],point=c(1,5:2,8:6,1),col=2)
#' }
#' ##3D example
#' \dontrun{
#' require(rgl)
#' data(nose)
#' points3d(shortnose.lm[1:9,])
#' lineplot(shortnose.lm[1:9,],point=list(c(1,3,2),c(3,4,5),c(8,6,5,7,9)),col=2)
#' }
#' 
#' @export
lineplot <- function(x,point,col=1,lwd=1,line_antialias = FALSE,lty=1,add=TRUE)
{
  
  if (dim(x)[2] == 3)
    {
      if (is.list(point)==TRUE)
        {
        start <- x[do.call(rbind, point)[, 1],]
        end<- x[do.call(rbind, point)[, 2], ]
        links <- do.call(rbind,Morpho::array2list(aperm(simplify2array(list(start,end)))))
        segments3d(links, lwd =2, col = 2)
        }
      else                                                
        {
          lines3d(x[point,1],x[point,2],x[point,3],col=col,lwd=lwd,line_antialias = line_antialias)
        }
    }
  else
    {
      if (!add)
        {
          plot(x,asp=1,cex=0,xlab="x-coordinate",ylab="y-coordinate")
        }
      if (is.list(point)==TRUE)
        {
          for (i in 1:length(point))
            {
              lines(x[point[[i]],1],x[point[[i]],2],col=col,lwd=lwd,lty=lty)
            }
        }
      else                                                
        {
          lines(x[point,1],x[point,2],col=col,lwd=lwd,lty=lty)
      }
    }
}
