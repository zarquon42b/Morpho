\name{meshDist}
\alias{meshDist.mesh3d}
\alias{meshDist.matrix}
\alias{meshDist}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  calculates and visualises distances between surface meshes or 3D
  coordinates and a surface mesh.
}
\description{
calculates the distances from a mesh or a set of 3D coordinates to another at each vertex; either
closest point or along the normals}
\usage{
\method{meshDist}{mesh3d}(x, mesh2=NULL,distvec=NULL, from = NULL, to
= NULL, steps = 20,ceiling = FALSE, file = "default", imagedim =
"100x800", uprange = 1, ray = FALSE, raytol = 50, save = FALSE, plot =
TRUE,sign=TRUE,tol=NULL,displace=FALSE,shade=TRUE,method=c("vcglib","morpho"),...)

\method{meshDist}{matrix}(x,mesh2=NULL,distvec=NULL,from=NULL,to=NULL,steps=20,ceiling=FALSE,uprange=1,plot=TRUE,sign=TRUE,tol=NULL,type=c("s","p"),radius=NULL,displace = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{reference mesh; object of class "mesh3d" or a n x 3 matrix
    containing 3D coordinates.
%%     ~~Describe \code{mesh1} here~~
}
  \item{mesh2}{target mesh: either object of class "mesh3d" or a
    character pointing to a surface mesh (ply, obj or stl file) 
%%     ~~Describe \code{mesh2} here~~
}
 \item{distvec}{vector: optional, a vector containing distances for each
   vertex of mesh1, if distvec != NULL, x will be ignored. 
%%     ~~Describe \code{mesh2} here~~
}
  \item{from}{numeric: minimum distance to be colorised; default is set to 0 mm
%%     ~~Describe \code{from} here~~
}
  \item{to}{numeric: maximum distance to be colorised; default is set to the maximum distance
%%     ~~Describe \code{to} here~~
}
  \item{steps}{integer: determines break points for color ramp: n steps
    will produce n-1 colors.
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
\item{sign}{logical: request signed distances. Only meaningful, if mesh2
  is specified or distvec contains signed distances.
%%     ~~Describe \code{raytol} here~~
}
\item{tol}{numeric: threshold to color distances within this threshold green.
}
\item{displace}{logical: if TRUE, displacement vectors between original
  and closest points are drawn colored according to the distance.
}
\item{shade}{logical: if FALSE, the rendering of the colored surface
  will be supressed.
}
\item{method}{accepts: "vcglib" and "morpho". vcglib uses a command line
  tool using vcglib headers, morpho uses fortran routines using kd-tree
  search for closest triangles.
}
\item{type}{character: "s" shows coordinates as spheres, while "p" shows
  3D dots.
}
\item{radius}{determines size of spheres; if not specified, optimal
  radius size will be estimated by centroid size of the configuration.
}
\item{\dots}{additional arguments passed to  \code{\link{shade3d}}. See
  \code{\link{rgl.material}} for details.
}
}
\details{this function needs the command line tools from the Auxiliaries section in \url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries} installed. 
}
\value{Returns an object of class "meshDist" if the input is a surface
  mesh and one of class "matrixDist" if input is a matrix containing 3D coordinates.
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{render.meshDist}},, \code{\link{export.meshDist}}, \code{\link{shade3d}}
}
\examples{
data(nose)##load data
##warp a mesh onto another landmark configuration:
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)


meshDist(warpnose.long,shortnose.mesh,method="m")



#use signed distances and
#color distances < 0.01 green:
meshDist(warpnose.long,shortnose.mesh,sign=TRUE,tol=0.01,method="m")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
