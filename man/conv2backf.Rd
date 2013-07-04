\name{conv2backf}
\alias{conv2backf}
\title{
invert faces' orientation of triangular mesh }
\description{
inverts faces' orientation of triangular mesh and recomputes vertex normals}
\usage{
conv2backf(mesh)
}
\arguments{
  \item{mesh}{
}
}

\value{returns resulting mesh
}

\author{Stefan Schlager
}


\seealso{
\code{\link{adnormals}}
}
\examples{
data(nose)
shade3d(shortnose.mesh,col=3)
noseinvert <- conv2backf(shortnose.mesh)
## show normals
plotNormals(noseinvert,long=0.01)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
