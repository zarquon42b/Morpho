% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeArea.r
\name{computeArea}
\alias{computeArea}
\title{Compute area enclosed within an irregular polygon}
\usage{
computeArea(x)
}
\arguments{
\item{x}{k x 2 or k x 3 matrix containing ordered coordinates forming the boundary of the area. For 3D-cases, the area should be closed to a 2D surface (see details below).}
}
\value{
returns a list containing
\item{area}{size of the enclosed area}
 \item{xpro2D}{projected coordinates of x in the 2D plane.}
 \item{poly}{object of class \code{sp} as defined by the \code{sp} package.}
\item{xpro3D}{For 3D-cases, this contains the projected coordinates of x rotated back into the original coordinate system}
}
\description{
Compute area enclosed within an irregular polygon - i.e. defined by curves
}
\details{
For 3D coordinates, a PCA is computed and only the first two PCs are used to compute the area. This is a projection of the coordinates onto a 2D plane spanned by those PCs.
}
\note{
in case custom planes are preferred, the data can first be projected onto such a custom defined plane via \code{\link{points2plane}} first.
}
\examples{
require(shapes)
require(sf)
myarea <- computeArea(gorf.dat[c(1,6:8,2:5),,1])
myarea$area
plot(myarea$poly)


## 3D example
data(boneData)
myarea3D <- computeArea(boneLM[c(4,2,3,7,5,6,8),,1])
plot(myarea3D$poly)
cent <- colMeans(myarea3D$xpro2D)
text(cent[1],cent[2],labels=paste0("Area=",round(myarea3D$area,digits=2)))
}
