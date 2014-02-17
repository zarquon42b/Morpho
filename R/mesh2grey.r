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
    


#' convert a colored mesh to greyscale.
#' 
#' convert the colors of a colored mesh to greyscale values
#' 
#' 
#' @param mesh Object of class mesh3d
#' @return returns a mesh with material$color replaced by greyscale rgb values.
#' @author Stefan Schlager
#' @seealso \code{\link{ply2mesh}},\code{\link{file2mesh}}
#' 
#' @export
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
        
