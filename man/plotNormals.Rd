\name{plotNormals}
\alias{plotNormals}
\title{plots the normals of a triangular surface mesh.}
\usage{
plotNormals(x, long = 1, lwd = 1, col = 1)
}
\arguments{
  \item{x}{object of class "mesh3d"}

  \item{long}{length of the normals (default is 1)}

  \item{lwd}{width of the normals}

  \item{col}{color of the normals}
}
\description{
visualises the vertex normals of a triangular surface mesh
of class mesh3d. If no normals are contained, they are
computed.
}
\examples{
\dontrun{
require(rgl)
data(nose)
plotNormals(shortnose.mesh,col=4,long=0.01)
shade3d(shortnose.mesh,col=3)
}
}
\author{
Stefan Schlager
}

