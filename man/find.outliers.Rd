\name{find.outliers}
\alias{find.outliers}
\title{Graphical interface to find outliers and/or to switch mislabeld landmarks 

}
\description{Graphical interface to find outliers and/or to switch mislabeld landmarks
}
\usage{
find.outliers(A, color = 4, lwd = 1, lcol = 2, mahalanobis=FALSE, PCuse=NULL)
}
\arguments{
  \item{A}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
}
  \item{color}{color of Landmarks points to be plotted
}
  \item{lwd}{linewidth visualizing distances of the individual landmarks from mean.
}
  \item{lcol}{color of lines visualizing distances of the individual landmarks from mean.
} 
  \item{mahalanobis}{logical: use mahalanobis distance to find outliers.
}  
\item{PCuse}{integer: Restrict mahalanobis distance to the first n
  Principal components.
} 
}
\details{This function performs a procrustes fit and tries sorts all
  specimen according to their distances (either Procrustes or
  Mahalanobis-distance) and provides visual help for rearranging
  landmarks and/or excluding outliers.

}
\value{

 \item{data.cleaned}{array (in original coordinate system) containing the changes applied and outliers eliminated }
 \item{outlier }{vector with integers indicating the positions in the
 original array that have been marked as outliers.}
 \item{dist.sort }{table showing the distance to mean for each observation - decreasing by distance}
 \item{type }{what kind of distance was used}
}


\author{Stefan Schlager
}


\seealso{
\code{\link{typprob}},\code{\link{typprobClass}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
