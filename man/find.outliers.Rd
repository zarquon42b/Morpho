\name{find.outliers}
\alias{find.outliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical interface to find outliers and/or to switch mislabeld landmarks 

}
\description{Graphical interface to find outliers and/or to switch mislabeld landmarks
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.outliers(A, color = 4, lwd = 1, lcol = 2, mahalanobis=FALSE, PCuse=NULL)
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
  \item{mahalanobis}{logical: use mahalanobis distance to find outliers.
%%     ~~Describe \code{lcol} here~~
}  
\item{PCuse}{integer: Restrict mahalanobis distance to the first n
  Principal components.
%%     ~~Describe \code{lcol} here~~
} 
}
\details{This function performs a procrustes fit and tries sorts all
  specimen according to their distances (either Procrustes or
  Mahalanobis-distance) and provides visual help for rearranging
  landmarks and/or excluding outliers.

}
\value{

%%  If it is a LIST, use
 \item{data.cleaned}{array (in original coordinate system) containing the changes applied and outliers eliminated }
 \item{outlier }{vector with integers indicating the positions in the
 original array that have been marked as outliers.}
 \item{dist.sort }{table showing the distance to mean for each observation - decreasing by distance}
 \item{type }{what kind of distance was used}
}


\author{Stefan Schlager
%%  ~~who you are~~
}


\seealso{
\code{\link{typprob}},\code{\link{typprobClass}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
