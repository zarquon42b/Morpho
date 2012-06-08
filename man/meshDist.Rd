\name{meshDist}
\alias{meshDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  calculates distances between surface meshes and renders and/or saves a
  colored mesh and image with color legend.
}
\description{
calculates the distances from a mesh to another at each vertex; either
closest point or along the normals}
\usage{
meshDist(mesh1, mesh2, from = NULL, to = NULL, steps = 20, ceiling = FALSE, file = "default", imagedim = "100x800", uprange = 1, ray = FALSE, raytol = 50, save = FALSE, plot = TRUE,sign=FALSE,tol=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh1}{reference mesh; object of class "mesh3d"
%%     ~~Describe \code{mesh1} here~~
}
  \item{mesh2}{target mesh: either object of class "mesh3d" or a
    character pointing to a surface mesh (ply, obj or stl file) 
%%     ~~Describe \code{mesh2} here~~
}
  \item{from}{numeric: minimum distance to color; default is set to 0 mm
%%     ~~Describe \code{from} here~~
}
  \item{to}{numeric: maximum distance to color; default is set to the maximum distance
%%     ~~Describe \code{to} here~~
}
  \item{steps}{integer: determines how many intermediate colors the color ramp has.
%%     ~~Describe \code{steps} here~~
}
  \item{ceiling}{logical: if TRUE, the next larger integer of "to" is
    used
%%     ~~Describe \code{ceiling} here~~
}
  \item{file}{character: filename for mesh and image files
    produced. E.g. "mydist" will produce the files mydist.ply and mydist.png
%%     ~~Describe \code{file} here~~
}
  \item{imagedim}{character of type 100x200 where 100 determines the
    width and 200 the height of the image.
%%     ~~Describe \code{imagedim} here~~
}
  \item{uprange}{numeric between 0 and 1: restricts "to" to a quantile of
    "to", if to is NULL.
%%     ~~Describe \code{uprange} here~~
}
  \item{ray}{logical: if TRUE, the search is along vertex normals.
%%     ~~Describe \code{ray} here~~
}
  \item{raytol}{maximum distance to follow a normal.
%%     ~~Describe \code{raytol} here~~
}
 \item{save}{logical: save a colored mesh.
%%     ~~Describe \code{raytol} here~~
}
\item{plot}{logical: visualise result as 3D-plot and distance charts
%%     ~~Describe \code{raytol} here~~
}
\item{sign}{logical: request signed distances.
%%     ~~Describe \code{raytol} here~~
}
\item{tol}{numeric: threshold to color distances within this threshold green.
}
}
\details{this function needs the command line tools from the Auxiliaries section in http://sourceforge.net/projects/morpho-rpackage/files/ installed. 
}
\value{
\item{colMesh}{object of mesh3d with colors added
}
\item{dists}{vector with distances.
}
\item{cols}{vector with color values.
}
\item{params}{list of parameters used.
  }
%%     ~~Describe \code{raytol} here~~
}
\references{
  Detection of inside/outside uses the algorithm proposed in:

  Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance Fields From Triangle Meshes. Informatics and Mathematical Modelling, .
}
\author{
Stefan Schlager}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{render.meshDist}}, \code{\link{export.meshDist}}
}
\examples{
data(nose)##load data
##warp a mesh onto another landmark configuration:
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)

## Not run:
meshDist(warpnose.long,shortnose.mesh)
## End(Not run)

## Not run:
#use signed distances and
#color distances < 0.01 green:
meshDist(warpnose.long,shortnose.mesh,sign=T,tol=0.01)
## End(Not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
