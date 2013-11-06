\name{write.pts}
\alias{write.pts}
\title{exports a matrix containing landmarks into .pts format}
\usage{
  write.pts(x, filename = dataname)
}
\arguments{
  \item{x}{k x m matrix containing landmark configuration}

  \item{filename}{character: Path/name of the requested
  output - extension will be added atuomatically. If not
  specified, the file will be named as the exported
  object.}
}
\description{
  exports a matrix containing landmarks into .pts format
  that can be read by IDAV Landmark.
}
\details{
  you can import the information into the program landmarks
  available at
  http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\examples{
data(nose)
write.pts(shortnose.lm, filename="shortnose")
}
\author{
  Stefan Schlager
}
\seealso{
  \code{\link{read.pts}}
}
\keyword{~kwd1}
\keyword{~kwd2}

