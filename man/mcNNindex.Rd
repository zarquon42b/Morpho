\name{mcNNindex}
\alias{mcNNindex}
\title{find nearest neighbours for point clouds}
\usage{
mcNNindex(target, query, cores = detectCores(), k = k, ...)
}
\arguments{
  \item{target}{\code{k x m} matrix containing data which
  to search.}

  \item{query}{\code{l x m} matrix containing data for
  which to search.}

  \item{cores}{integer: amount of CPU-cores to be used.
  Speed benefits are only relevant for \code{k > 20}}

  \item{k}{integer: how many closest points are sought.}

  \item{\dots}{additional arguments - currently unused.}
}
\value{
  \code{l x k } matrix containing indices of closest
  points.
}
\description{
  find nearest neighbours for point clouds by using
  algorithms from the ANN library. This is just a wrapper
  for the function ann from the package yaImpute, enabling
  parallel processing.
}
\details{
  wraps the function \code{ann} from package 'yaImpute' to
  allow multicore processing
}
\examples{
require(rgl)
data(nose)
# find closest vertex on surface for each landmark
clost <- mcNNindex(vert2points(shortnose.mesh),shortnose.lm, k=1,
mc.cores=1)
\dontrun{
spheres3d(vert2points(shortnose.mesh)[clost,],col=2,radius=0.3)
spheres3d(shortnose.lm,radius=0.3)
wire3d(shortnose.mesh)
}
}
\seealso{
  \code{\link{closemeshKD}}
}
\keyword{~kwd1}
\keyword{~kwd2}

