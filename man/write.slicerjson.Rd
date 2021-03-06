% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fiducials.r
\name{write.slicerjson}
\alias{write.slicerjson}
\title{Export landmarks (or any 3D coordinates) to the new slicer json format}
\usage{
write.slicerjson(
  x,
  filename = dataname,
  type = c("Fiducial", "Curve", "ClosedCurve"),
  coordinateSystem = c("LPS", "RAS"),
  labels = dataname
)
}
\arguments{
\item{x}{k x 3 matrix containing 3D coordinates}

\item{filename}{will be substituted with ".mrk.json"}

\item{type}{character: specify type of coordinates. Can be any of "Fiducial", "Curve", "ClosedCurve".}

\item{coordinateSystem}{character: specify coordinate system the data are in. Can be "LPS" or "RAS".}

\item{labels}{character or character vector containing landmark labels.}
}
\description{
Export landmarks (or any 3D coordinates) to the new slicer json format
}
