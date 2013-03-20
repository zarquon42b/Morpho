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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{sample data: matrix or vector
%%     ~~Describe \code{x} here~~
}
  \item{output}{logical: if TRUE results are returned
%%     ~~Describe \code{output} here~~
}
  \item{square}{plot in a square window - outliers might be cut off.
%%     ~~Describe \code{square} here~~
}
}

\value{if output=TRUE, the following values are returned
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{x }{Distances from an expected Gaussian distribution}
  \item{y }{Observed distances - sorted}
   \item{d }{Observed distances - unsorted}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Stefan Schlager
}


\seealso{
 \code{\link{qqplot}}
}
\examples{
### create normally distributed data
data <- mvrnorm(100,mu=rep(0,5),Sigma = diag(5:1))
qqmat(data)

###create non normally distributed data
data1 <- rchisq(100,df=3)
qqmat(data1,square=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
