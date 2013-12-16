\name{updateNormals}
\alias{facenormals}
\alias{updateNormals}
\title{Compute face or vertex normals of a triangular mesh}
\usage{
updateNormals(x, angle = TRUE)

facenormals(x)
}
\arguments{
  \item{x}{triangular mesh of class "mesh3d"}

  \item{angle}{logical: if TRUE, angle weighted normals are
  used.}
}
\value{
adnormals returns mesh with updated vertex normals.

facenormals returns an object of class "mesh3d". With
\item{vb }{faces' barycenters} \item{normals }{faces'
normals}
}
\description{
Compute face or vertex normals of a triangular mesh of
class "mesh3d"
}
\note{
only supports triangular meshes
}
\examples{
require(rgl)
require(Morpho)
data(nose)
### calculate vertex normals
shortnose.mesh$normals <- NULL ##remove normals
\dontrun{
shade3d(shortnose.mesh,col=3)##render
}
shortnose.mesh <- updateNormals(shortnose.mesh)
\dontrun{
rgl.clear()
shade3d(shortnose.mesh,col=3)##smoothly rendered now
}
## calculate facenormals
facemesh <- facenormals(shortnose.mesh)
plotNormals(facemesh,long=0.01)
\dontrun{
points3d(vert2points(facemesh),col=2)
wire3d(shortnose.mesh)
}
}
\author{
Stefan Schlager
}
\references{
Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating
Signed Distance Fields From Triangle Meshes. Informatics
and Mathematical Modelling, .
}
\seealso{
\code{\link{ply2mesh}}
}
\keyword{~kwd1}
\keyword{~kwd2}

