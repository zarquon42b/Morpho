\name{kendalldist}
\alias{kendalldist}
\title{Calculates the Riemannian distance between two superimposed
  landmark configs.
}
\description{Calculates the Riemannian distance between two superimposed
  landmark configs.
}
\usage{
kendalldist(x, y)
}
\arguments{
  \item{x}{Matrix containing landmark coordinates.
  }
  \item{y}{Matrix containing landmark coordinates.
  }
}
\value{returns Riemannian distance
}
\examples{
library(shapes)
OPA <- rotonto(gorf.dat[,,1],gorf.dat[,,2])
kendalldist(OPA$X,OPA$Y)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
