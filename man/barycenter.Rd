\name{barycenter}
\alias{barycenter}
\title{
  calculates the barycenters for all faces of a triangular mesh
}
\description{
  calculates the barycenters for all faces of a triangular mesh
}
\usage{
barycenter(mesh)
}

\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'
  }
}

\value{k x 3 matrix of barycenters for all \code{k} faces of input mesh.
}

\seealso{
  \code{\link{closemeshKD}}
  }
\examples{
data(nose)
bary <- barycenter(shortnose.mesh)
##visualize mesh
wire3d(shortnose.mesh)
# visualize barycenters
points3d(bary, col=2)
## now each triangle is equipped with a point in its barycenter
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
