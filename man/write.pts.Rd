\name{write.pts}
\alias{write.pts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
exports a matrix containing landmarks into .pts format
}
\description{exports a matrix containing landmarks into .pts format that can be read by IDAV Landmark.

}
\usage{
write.pts(X, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{k x m matrix containing landmark configuration
%%     ~~Describe \code{X} here~~
}
  \item{filename}{ filename - extension .pts will be added automatically
%%     ~~Describe \code{filename} here~~
}
}
\details{you can import the information into the program landmarks available at http://graphics.idav.ucdavis.edu/research/EvoMorph
}

\author{  Stefan Schlager
%%  ~~who you are~~
}


\seealso{
\code{\link{read.pts}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
