RGB2Grey <- function(x,coefs = c(0.3, 0.59, 0.11))
  {
    matr <- FALSE
    if (is.matrix(x))
      {dimx <- dim(x)
       x <- as.vector(x)
       matr <- TRUE
     }
    rgbval <- col2rgb(x)
    greyval <- round((t(coefs)%*%rgbval))
    greyval <- rgb(greyval,greyval,greyval,maxColorValue=255)

    if (matr)
      {
        greyval <- matrix(greyval,dimx[1],dimx[2])
      }
    return(greyval)
  }
    
mesh2grey <- function(mesh)
  {
    if (is.null(mesh$material$color))
      {
        stop("mesh contains no colors")
      }
    else
      {
        mesh$material$color <- RGB2Grey(mesh$material$color)
      }
    return(mesh)
  }
        
