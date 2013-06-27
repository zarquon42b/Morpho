\name{scalemesh}
\alias{scalemesh}
\title{scale a mesh of class "mesh3d"

}
\description{scales (the vertices of a mesh by a scalar
}
\usage{
scalemesh(mesh, size, center = c("bbox", "mean", "none"))
}
\arguments{
  \item{mesh}{object of class "mesh3d"
}
  \item{size}{numeric: scale factor
}
\item{center}{character: method to position center of mesh after
  scaling: values are "bbox", and "mean". See Details for more info.
}
}
\details{The mesh's center is determined either as mean of the bounding
  box (center="bbox") or mean of vertex coordinates (center="mean") and
  then scaled according to the scaling factor. If center="none", vertex
  coordinates will simply be multiplied by "size".
}
\value{returns a scaled mesh

}

\author{
Stefan Schlager}


\seealso{
\code{\link{rotmesh.onto}}
}
\examples{
data(nose)
#inflate mesh by factor 4
largenose <- scalemesh(shortnose.mesh,4)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
