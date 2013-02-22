\name{rotaxis3d}
\alias{rotaxis3d}
\alias{rotaxis3d.matrix}
\alias{rotaxis3d.mesh3d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rotate an object around an arbitrary axis in 3D}
\description{
Rotate an object (matrix or triangular mesh) around an 3D-axis defined
by two points.
}
\usage{
\method{rotaxis}{matrix}(x, pt1, pt2 = c(0, 0, 0), theta)
\method{rotaxis}{mesh3d}(x, pt1, pt2 = c(0, 0, 0), theta)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{k x 3 matrix containing 3D-coordinates or a triangular mesh
    of class "mesh3d".
}
  \item{pt1}{numeric vector of length , defining first point on the
    rotation axis.
}
  \item{pt2}{numeric vector of length , defining second point on the
    rotation axis.
}
  \item{theta}{angle to rotate in radians.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  returns rotated object (including updated normals for mesh3d objects)

}
\references{
  http://en.wikipedia.org/wiki/Rotation_matrix
}
\author{
Stefan Schlager}



\seealso{
 \code{\link{rotonto}}, \code{\link{rotmesh.onto}}
}
\examples{
data(nose)
shade3d(shortnose.mesh,col=3,specular=1)
shrot.rot <- rotaxis3d(shortnose.mesh,pt1=c(1,1,1),theta=pi)
shade3d(shrot.rot,col=2)
###print rotation axis

lines3d(rbind(rep(-0.1,3),rep(0.1,3)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
