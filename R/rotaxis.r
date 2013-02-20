rotaxis <- function(u,theta)
  {
    crossmat <- function(x)
      {out <- matrix(c(0,x[3],-x[2],-x[3],0,x[1],x[2],-x[1],0),3,3)
     }
       
    u <- u/sqrt(sum(crossprod(u)))
    I <- diag(rep(1,3))
    R <- I*cos(theta)+sin(theta)*crossmat(u)+(1-(cos(theta)))*tcrossprod(u)
    return(R)
  }

rotaxisPoint <- function(mat,x,y,theta)
  {
    u <- y-x
    ### translate axis
    matrixTrans <- t(t(mat)-x)
    rotmatrix <- rotaxis(u,theta)
    out <- t(t(matrixTrans%*%rotmatrix)+x)
    return(out)
  }

rotaxisMesh <- function(mesh,x,y,theta)
  {
    mat <- vert2points(mesh)
    vb <- rotaxisPoint(mat,x,y,theta)
    mesh$vb[1:3,] <- t(vb)
    mesh <- adnormals(mesh)
    invisible(mesh)
  }
