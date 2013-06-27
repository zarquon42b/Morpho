\name{rotonmat}
\alias{rotonmat}
\title{
rotate matrix of landmarks}
\description{
rotate matrix of landmarks by using a rotation determined by two matrices.
}
\usage{
rotonmat(X, refmat, tarmat, scale = TRUE, reflection = FALSE)
}
\arguments{
  \item{X}{Matrix to be rotated
}
  \item{refmat}{reference matrix used to estimate rotation parameters
}
  \item{tarmat}{target matrix used to estimate rotation parameters
}
  \item{scale}{logical: requests scaling to minimize sums of squared distances
}
  \item{reflection}{logical: if TRUE, reflections are allowed.}
}
\details{A matrix is rotated by rotation parameters determined by two
  different matrices. This is usefull, if the rotation parameters are to
  be estimated by a subset of landmark coordinates.
}
\value{returns rotated X.

}

\author{
Stefan Schlager}


\seealso{
\code{\link{rotonto}},\code{\link{rotmesh.onto}}
}
\examples{

data(nose)
shortnose.rot <-
rotonmat(shortnose.lm,shortnose.lm[1:9,],longnose.lm[1:9,])

##view result
deformGrid3d(shortnose.rot,shortnose.lm,ngrid=0)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
