% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fiducials.r
\name{read.fcsv}
\alias{read.fcsv}
\title{read fiducials from slicer4}
\usage{
read.fcsv(x, na = NULL, lps2ras = FALSE)
}
\arguments{
\item{x}{filename}

\item{na}{value to be replaced by NA}

\item{lps2ras}{logical: if the coordinate system is LPS and \code{lps2ras=TRUE}, the data will be rotated into the RAS space by inverting the first two dimensions using \code{\link{LPS2RAS}}.}
}
\value{
a k x 3 matrix with landmarks
}
\description{
read fiducials from slicer4
}
