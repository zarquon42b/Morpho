\name{qqmat}
\alias{qqmat}
\title{Q-Q plot to assess normality of data
}
\description{
qqmat plots Mahalanobisdistances of a given sample to those expected
from a Gaussian distribution}
\usage{
qqmat(x, output = FALSE, square = TRUE)
}
\arguments{
  \item{x}{sample data: matrix or vector
}
  \item{output}{logical: if TRUE results are returned
}
  \item{square}{plot in a square window - outliers might be cut off.
}
}

\value{if output=TRUE, the following values are returned
  \item{x }{Distances from an expected Gaussian distribution}
  \item{y }{Observed distances - sorted}
   \item{d }{Observed distances - unsorted}
}

\author{
  Stefan Schlager
}


\seealso{
 \code{\link{qqplot}}
}
\examples{
require(MASS)
### create normally distributed data
data <- mvrnorm(100,mu=rep(0,5),Sigma = diag(5:1))
qqmat(data)

###create non normally distributed data
data1 <- rchisq(100,df=3)
qqmat(data1,square=FALSE)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
