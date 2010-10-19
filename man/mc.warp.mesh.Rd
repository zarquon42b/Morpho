\name{mc.warp.mesh}
\alias{mc.warp.mesh}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Multicore Version of warp.mesh
}
\description{
      warps an the surface of a mesh3d object onto another configuration via reference and target landmark configuration.
}
\usage{
mc.warp.mesh(mesh, matr, matt, updateNormals = TRUE) 

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
    \code{\link{ply2mesh},\link{file2mesh},\link{mesh2ply}}
}
