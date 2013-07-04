\name{write.pts}
\alias{write.pts}
\title{
exports a matrix containing landmarks into .pts format
}
\description{exports a matrix containing landmarks into .pts format that can be read by IDAV Landmark.

}
\usage{
write.pts(X, filename)
}
\arguments{
  \item{X}{k x m matrix containing landmark configuration
}
  \item{filename}{ filename - extension .pts will be added automatically
}
}
\details{you can import the information into the program landmarks available at http://graphics.idav.ucdavis.edu/research/EvoMorph
}

\author{  Stefan Schlager
}


\seealso{
\code{\link{read.pts}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
