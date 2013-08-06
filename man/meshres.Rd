\name{meshres}
\alias{meshres}
\alias{meshres.mesh3d}
\title{calculate average edge length of a triangular mesh}
\usage{
  \method{meshres}{mesh3d}(mesh)
}
\arguments{
  \item{mesh}{triangular mesh stored as object of class
  "mesh3d"}
}
\value{
 returns average edge length (a.k.a. mesh resolution)
 }
\description{
  calculate average edge length of a triangular mesh, by
  iterating over all faces.
}
\examples{
data(boneData)
mres <- meshres(skull_0144_ch_fe.mesh)

}
\author{
  Stefan Schlager
}
\keyword{~kwd1}
\keyword{~kwd2}

