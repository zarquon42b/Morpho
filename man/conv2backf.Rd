\name{conv2backf}
\alias{conv2backf}
\title{invert faces' orientation of triangular mesh}
\usage{
conv2backf(mesh)
}
\arguments{
  \item{mesh}{triangular mesh of class \code{mesh3d}}
}
\value{
returns resulting mesh
}
\description{
inverts faces' orientation of triangular mesh and
recomputes vertex normals
}
\examples{
require(rgl)
data(nose)
\dontrun{
shade3d(shortnose.mesh,col=3)
}
noseinvert <- conv2backf(shortnose.mesh)
## show normals
plotNormals(noseinvert,long=0.01)
}
\author{
Stefan Schlager
}
\seealso{
\code{\link{adnormals}}
}
\keyword{~kwd1}
\keyword{~kwd2}

