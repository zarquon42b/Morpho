\name{vecx}
\alias{vecx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{convert an 3D array into a matrix
}
\description{converts a 3D-array (e.g. containing landmark coordinates)
  into a matrix, one row per specimen.
}
\usage{
vecx(x, byrow = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{array
    %%     ~~Describe \code{x} here~~
  }
  \item{byrow}{logical: if TRUE, the resulting vector for each specimen
    will be x1,y1,z1,x2,y2,z2,..., and
    x1,x2,...,y1,y2,...,z1,z2,... otherwise (default)
    %%     ~~Describe \code{byrow} here~~
  }
}


\value{returns a matrix with one row per specimen
}

\author{
  Stefan Schlager}

\examples{
data <- vecx(gorf.dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
