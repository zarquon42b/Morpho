\name{PCdist}
\alias{PCdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 correlation between a reduced space and the original space
}
\description{
Calculates the correlation between distances in a reduced space and the original space
}
\usage{
PCdist(PCs, PCscores, x = 5, plot.type="b")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PCs}{m x k matrix of Principal Components where m is the k is the number of PCs.
%%     ~~Describe \code{PCs} here~~
}
  \item{PCscores}{n x m matrix of Principal Component scores where n is the number of observations.
%%     ~~Describe \code{PCscores} here~~
}
  \item{x}{integer: increment for every x-th PC the subspace to fullspace correlation will be calculated.
%%     ~~Describe \code{x} here~~
}
 \item{plot.type }{"b"=barplot of correlation values, "s"=line between correlation values.
%%     ~~Describe \code{x} here~~
}
}

\value{a vector of R-squared values between subspace and fullspace distances and a barplot depicting the correlations belonging to the subspace.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Stefan Schlager}

\examples{
library(shapes)
a <- procSym(gorf.dat)
PCdist(a$PCs, a$PCscores, x = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
