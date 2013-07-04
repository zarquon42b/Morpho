\name{vecx}
\alias{vecx}
\title{convert an 3D array into a matrix
}
\description{converts a 3D-array (e.g. containing landmark coordinates)
  into a matrix, one row per specimen.
}
\usage{
vecx(x, byrow = FALSE)
}
\arguments{
  \item{x}{array
    
  }
  \item{byrow}{logical: if TRUE, the resulting vector for each specimen
    will be x1,y1,z1,x2,y2,z2,..., and
    x1,x2,...,y1,y2,...,z1,z2,... otherwise (default)
    
  }
}


\value{returns a matrix with one row per specimen
}

\author{
  Stefan Schlager}

\examples{
library(shapes)
data <- vecx(gorf.dat)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
