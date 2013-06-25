\name{scalemesh}
\alias{scalemesh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{scale a mesh of class "mesh3d"

}
\description{scales (the vertices of a mesh by a scalar
}
\usage{
scalemesh(mesh, size, center = c("bbox", "mean", "none"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{object of class "mesh3d"
%%     ~~Describe \code{mesh} here~~
}
  \item{size}{numeric: scale factor
%%     ~~Describe \code{size} here~~
}
\item{center}{character: method to position center of mesh after
  scaling: values are "bbox", and "mean". See Details for more info.
%%     ~~Describe \code{size} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
