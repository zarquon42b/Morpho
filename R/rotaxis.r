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

rotaxis3d <- function(x,pt1,pt2=c(0,0,0),theta) UseMethod("rotaxis3d")

rotaxis3d.matrix <- function(x,pt1,pt2=c(0,0,0),theta)
  {
    u <- pt2-pt1
    ### translate axis
    matrixTrans <- t(t(x)-pt1)
    rotmatrix <- rotaxisMat(u,theta)
    out <- t(t(matrixTrans%*%rotmatrix)+pt1)
    return(out)
  }

rotaxis3d.mesh3d <- function(x,pt1,pt2=c(0,0,0),theta)
  {
    mat <- vert2points(x)
    vb <- rotaxis3d(mat,pt1,pt2,theta)
    x$vb[1:3,] <- t(vb)
    x <- adnormals(x)
    invisible(x)
  }
