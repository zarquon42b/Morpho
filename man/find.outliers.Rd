\name{find.outliers}
\alias{find.outliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical interface to find outliers and/or to switch mislabeld landmarks 

}
\description{Graphical interface to find outliers and/or to switch mislabeld landmarks
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.outliers(A, color = 4, lwd = 1, lcol = 2, sz = 0, output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
%%     ~~Describe \code{A} here~~
}
  \item{color}{color of Landmarks points to be plotted
%%     ~~Describe \code{color} here~~
}
  \item{lwd}{linewidth visualizing distances of the individual landmarks from mean.
%%     ~~Describe \code{lwd} here~~
}
  \item{lcol}{color of lines visualizing distances of the individual landmarks from mean.
%%     ~~Describe \code{lcol} here~~
}

}
\details{

}
\value{

%%  If it is a LIST, use
 \item{data.cleaned}{array (in original coordinate system) containing the changes applied and outliers eliminated }
 \item{outlier }{Description of 'comp2'}
 \item{dist.sort }{table showing the distance to mean for each observation - decreasing by distance}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Stefan Schlager
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
