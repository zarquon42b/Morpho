\name{deformGrid3d}
\alias{deformGrid3d}
\title{
visualise differences between two superimposed sets of 3D landmarks}
\description{visualise differences between two superimposed sets of 3D
  landmarks by deforming a cubic grid based on a thin-plate spline interpolation}
  
\usage{
deformGrid3d(matrix, tarmatrix, ngrid = 10, lwd = 1, showaxis = c(1, 2,
3),both = T, lines = TRUE, lcol = 1, add = FALSE, col1 = 2,
col2 = 3, type = c("s", "p"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{reference matrix containing 3D landmark coordinates.
}
  \item{tarmatrix}{target matrix containing 3D landmark coordinates.
%%     ~~Describe \code{tarmatrix} here~~
}
  \item{ngrid}{number of grid lines to be plotted; ngrid=0 suppresses
    grid creation.
%%     ~~Describe \code{ngrid} here~~
}
  \item{lwd}{width of lines connecting landmarks.
%%     ~~Describe \code{lwd} here~~
}
  \item{showaxis}{which dimensions of the grid to be plotted.
%%     ~~Describe \code{showaxis} here~~
}
  \item{both}{logical: if FALSE, only "matrix" will be plotted.
%%     ~~Describe \code{both} here~~
}
  \item{lines}{logical: if TRUE, lines between landmarks will be plotted.
%%     ~~Describe \code{lines} here~~
}
  \item{lcol}{color of lines
%%     ~~Describe \code{lcol} here~~
}
  \item{add}{logical: add to existing rgl window.
%%     ~~Describe \code{add} here~~
}
  \item{col1}{color of "matrix"
%%     ~~Describe \code{col1} here~~
}
  \item{col2}{color of "tarmat"
%%     ~~Describe \code{col2} here~~
}
  \item{type}{"s" renders landmarks as spheres; "p" as points - much
    faster for very large pointclouds.
%%     ~~Describe \code{type} here~~
}
}


\author{
Stefan Schlager}


\seealso{\code{\link{tps3d}}
}
\examples{
data(nose)
\dontrun{
deformGrid3d(shortnose.lm,longnose.lm,ngrid=10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
