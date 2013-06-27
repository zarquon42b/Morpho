\name{covW}
\alias{covW}
\title{calculate the pooled within groups covariance matrix

}
\description{calculate the pooled within groups covariance matrix
}
\usage{
covW(data, groups)
}
\arguments{
  \item{data}{a matrix containing data
}
  \item{groups}{grouping variables.

}
}

\value{
  Returns the pooled within group covariance matrix.
}
\author{
Stefan Schlager}

\seealso{
 \code{\link{cov}}, \code{\link{typprobClass}}
}
\examples{
data(iris)
poolCov <- covW(iris[,1:4],iris[,5])
}
