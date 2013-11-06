#' calculate a rotation matrix around an arbitrary axis through the origin in
#' 3D
#' 
#' calculate a rotation matrix around an arbitrary axis in 3D
#' 
#' @param u a vector around which to rotate
#' @param theta angle in radians to rotate
#' @return returns 3x3 rotation matrix
#' @seealso \code{\link{rotaxis3d}}
#' @references http://en.wikipedia.org/wiki/Rotation_matrix
#' @keywords ~kwd1 ~kwd2
#' @export rotaxisMat
rotaxisMat <- function(u,theta)
  {
    crossmat <- function(x)
      {
          out <- matrix(c(0,x[3],-x[2],-x[3],0,x[1],x[2],-x[1],0),3,3)
          return(out)
      }
       
    u <- u/sqrt(sum(crossprod(u)))
    I <- diag(rep(1,3))
    R <- I*cos(theta)+sin(theta)*crossmat(u)+(1-(cos(theta)))*tcrossprod(u)
    return(R)
  }



#' Rotate an object around an arbitrary axis in 3D
#' 
#' Rotate an object (matrix or triangular mesh) around an 3D-axis defined by
#' two points.
#' 
#' @title  Rotate an object (matrix or mesh) around an arbitrary axis in 3D
#' @param x k x 3 matrix containing 3D-coordinates or a triangular mesh of
#' class "mesh3d".
#' @param pt1 numeric vector of length 3, defining first point on the rotation
#' axis.
#' @param pt2 numeric vector of length 3, defining second point on the rotation
#' axis.
#' @param theta angle to rotate in radians. With pt1 being the viewpoint, the
#' rotation is counterclockwise.
#' @return returns rotated object (including updated normals for mesh3d
#' objects)
#' @author Stefan Schlager
#' @seealso \code{\link{rotonto}}, \code{\link{rotmesh.onto}}
#' @references http://en.wikipedia.org/wiki/Rotation_matrix
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' shrot.rot <- rotaxis3d(shortnose.mesh,pt1=c(1,1,1),theta=pi)
#' \dontrun{
#' shade3d(shortnose.mesh,col=3,specular=1)
#' shade3d(shrot.rot,col=2)
#' 
#' ###print rotation axis
#' #' lines3d(rbind(rep(-0.1,3),rep(0.1,3)))
#' }
#' @export rotaxis3d
rotaxis3d <- function(x,pt1,pt2=c(0,0,0),theta) UseMethod("rotaxis3d")

#' @rdname rotaxis3d
#' @method rotaxis3d matrix
#' @S3method rotaxis3d matrix
rotaxis3d.matrix <- function(x,pt1,pt2=c(0,0,0),theta)
  {
    u <- pt2-pt1
    ### translate axis
    matrixTrans <- t(t(x)-pt1)
    rotmatrix <- rotaxisMat(u,theta)
    out <- t(t(matrixTrans%*%rotmatrix)+pt1)
    return(out)
  }
#' @rdname rotaxis3d
#' @method rotaxis3d mesh3d
#' @S3method rotaxis3d mesh3d
rotaxis3d.mesh3d <- function(x,pt1,pt2=c(0,0,0),theta)
  {
    mat <- vert2points(x)
    vb <- rotaxis3d(mat,pt1,pt2,theta)
    x$vb[1:3,] <- t(vb)
    x <- adnormals(x)
    invisible(x)
  }
