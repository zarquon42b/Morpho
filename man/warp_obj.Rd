\name{warp_obj}
\Rdversion{1.1}
\alias{warp_obj}

\title{
warping a mesh of an imported obj-file onto another configuration
}
\description{
      warps an the surface of an imported obj-file onto another configuration via reference and target landmark configuration and allows to save it to an .obj file.
}
\usage{
warp_obj(obj, matr, matt,save=TRUE,filename="default")
}

\arguments{
  
  \item{obj}{object imported via "read.obj" command

  }
  \item{matr}{matrix of landmarks on the reference surface

  }
  \item{matt}{matrix of corresponding landmarks on the target surface

  }
  \item{save}{Logical: requests the output to be saved into an obj-file in the working directory
  }
  \item{filename}{character string to determine the name of the obj-file. The fileextension will be added automatically (filename="surface" will create a file named surface.obj)
  }
}





\details{
    the surface is mapped via the tps3d function onto the target shape.
}
\value{
  \item{obj}{dataframe in .obj-form, which can be used to render with the function render.obj.}
}


\author{
Stefan Schlager
}

\seealso{
    \code{\link{render.obj}}
}