\name{deformGrid3d}
\alias{deformGrid3d}
\title{visualise differences between two superimposed sets of 3D landmarks}
\usage{
deformGrid3d(matrix, tarmatrix, ngrid = 10, lwd = 1, showaxis = c(1, 2),
  both = T, lines = TRUE, lcol = 1, add = FALSE, col1 = 2, col2 = 3,
  type = c("s", "p"))
}
\arguments{
  \item{matrix}{reference matrix containing 3D landmark
  coordinates.}

  \item{tarmatrix}{target matrix containing 3D landmark
  coordinates.}

  \item{ngrid}{number of grid lines to be plotted; ngrid=0
  suppresses grid creation.}

  \item{lwd}{width of lines connecting landmarks.}

  \item{showaxis}{integer (vector): which dimensions of the
  grid to be plotted. Options are combinations of 1,2 and
  3.}

  \item{both}{logical: if FALSE, only "matrix" will be
  plotted.}

  \item{lines}{logical: if TRUE, lines between landmarks
  will be plotted.}

  \item{lcol}{color of lines}

  \item{add}{logical: add to existing rgl window.}

  \item{col1}{color of "matrix"}

  \item{col2}{color of "tarmat"}

  \item{type}{"s" renders landmarks as spheres; "p" as
  points - much faster for very large pointclouds.}
}
\description{
visualise differences between two superimposed sets of 3D
landmarks by deforming a cubic grid based on a thin-plate
spline interpolation
}
\examples{
\dontrun{
data(nose)
deformGrid3d(shortnose.lm,longnose.lm,ngrid=10)
}
}
\author{
Stefan Schlager
}
\seealso{
\code{\link{tps3d}}
}

