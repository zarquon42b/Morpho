\name{warp.mesh}
\alias{warp.mesh}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
warping a mesh onto another configuration
}
\description{
      warps an the surface of a mesh3d object onto another configuration
      via reference and target landmark configuration by using a
      thin-plate spline interpolation.
}
\usage{
warp.mesh(mesh, matr, matt, lambda = 0, updateNormals = TRUE) 

}

\arguments{
  
  \item{mesh}{object of class "mesh3d"

  }
  \item{matr}{matrix of landmarks on the reference surface

  }
  \item{matt}{matrix of corresponding landmarks on the target surface

  }
  \item{updateNormals}{Logical: requests the (re)calculation of vertex normals.
  }
  \item{lambda}{integer: regularisation parameter of the TPS.

}
}





\details{
    the surface is mapped via the tps3d function onto the target shape.
}
\value{
  \item{mesh}{list of class "mesh3d"}
}


\author{
Stefan Schlager
}

\seealso{
    \code{\link{ply2mesh},\link{file2mesh},\link{mesh2ply},\link{warpmovie3d}}
  }
  
  \examples{
data(nose)##load data
##warp a mesh onto another landmark configuration:
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
shade3d(warpnose.long,col=skin1)

}
