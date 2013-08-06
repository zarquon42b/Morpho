\name{mcNNindex}
\alias{mcNNindex}
\title{
find nearest neighbours for point clouds
}
\description{find nearest neighbours for point clouds by using
  algorithms from the ANN library. This is just a wrapper for the
  function ann from the package yaImpute, enabling parallel processing.

}
\usage{
mcNNindex(target, query, cores = detectCores(), k = k, ...)
}
\arguments{
  \item{target}{\code{k x m} matrix containing data which to search.
}
  \item{query}{\code{l x m} matrix containing data for which to search.
}
  \item{cores}{integer: amount of CPU-cores to be used. Speed benefits
    are only relevant for \code{k > 20}
}
  \item{k}{integer: how many closest points are sought.
}
  \item{\dots}{additional arguments - currently unused.

  }
}
\details{
wraps the function \code{ann} from package 'yaImpute' to allow multicore
processing}
\value{ \code{l x k } matrix containing indices of closest points.
}

\seealso{ \code{\link{closemeshKD}}
}
\examples{
data(boneData)
# find closest vertex on surface for each landmark
clost <- mcNNindex(vert2points(skull_0144_ch_fe.mesh),boneLM[,,1], k=1, mc.cores=1)
spheres3d(vert2points(skull_0144_ch_fe.mesh)[clost,],col=2,radius=0.3)
spheres3d(boneLM[,,1],radius=0.3)
wire3d(skull_0144_ch_fe.mesh)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
