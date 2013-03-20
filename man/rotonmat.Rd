\name{rotonmat}
\alias{rotonmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rotate matrix of landmarks}
\description{
rotate matrix of landmarks by using a rotation determined by two matrices.
}
\usage{
rotonmat(X, refmat, tarmat, scale = TRUE, reflection = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Matrix to be rotated
}
  \item{refmat}{reference matrix used to estimate rotation parameters
%%     ~~Describe \code{refmat} here~~
}
  \item{tarmat}{target matrix used to estimate rotation parameters
}
  \item{scale}{logical: requests scaling to minimize sums of squared distances
%%     ~~Describe \code{scale} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
