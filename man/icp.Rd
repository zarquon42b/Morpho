\name{icp}
\alias{icp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Iterative closest point matching between two triangular meshes.
%%  ~~function to do ... ~~
}
\description{
performs rigid body transformations (and scaling if requested) to map a
reference mesh onto a target mesh}
\usage{
icp(mesh1, mesh2, iterations = 3, scale = T, lm1 = NULL, lm2 = NULL, uprange = 1, rhotol = pi, k = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh1}{object of class "mesh3d".
%%     ~~Describe \code{mesh1} here~~
}
  \item{mesh2}{object of class "mesh3d".
%%     ~~Describe \code{mesh2} here~~
}
  \item{iterations}{integer
%%     ~~Describe \code{iterations} here~~
}
  \item{scale}{lgoical: scale reference mesh to optimize fit.
%%     ~~Describe \code{scale} here~~
}
  \item{lm1}{optional: kx3 matrix containing reference points on mesh1.
%%     ~~Describe \code{lm1} here~~
}
  \item{lm2}{optional: kx3 matrix containing reference points on mesh2.
%%     ~~Describe \code{lm2} here~~
}
  \item{uprange}{quantile of distances between vertices of mesh1 and
    closest points on mesh2. All hit points on mesh2 farther away than
    the specified quantile are not used for matching.
%%     ~~Describe \code{uprange} here~~
}
  \item{rhotol}{maximum allowed angle of which normals are allowed to
    differ between reference points and hit target points. Hit points
    above this threshold will be neglected.
%%     ~~Describe \code{rhotol} here~~
}
  \item{k}{integer: number of closest triangles on target that will be
    searched for hits.
%%     ~~Describe \code{k} here~~
}
}
\details{
the registration is done by minimising squared distances between
reference vertices and closest points on the target surface (a.k.a. Procrustes registration)}
\value{
Returns the rotated mesh1.
}
\references{Zhang Z. 1994. Iterative point matching for registration of free-form curves and surfaces International Journal of Computer Vision 13:119-152.
}
\author{
Stefan Schlager}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{rotmesh.onto}},  \code{\link{rotonto}}
}
\examples{
data(nose)
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
rotnose <- icp(warpnose.long,shortnose.mesh,lm1=longnose.lm,lm2=shortnose.lm,rhotol=0.7,uprange=0.9)
shade3d(rotnose,col=2,alpha=0.7)
shade3d(shortnose.mesh,col=3,alpha=0.7)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
