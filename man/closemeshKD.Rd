\name{closemeshKD}
\alias{closemeshKD}
\title{Project coordinates onto a target triangular surface mesh.
}

\description{For a set of 3D-coordinates the closest matches on a target
  surface are determined and normals at as well as distances to that point are calculated. 
}
\usage{
closemeshKD(x, mesh, k = 50, sign = FALSE, cores = 1, method = 0, ...)
}
\arguments{
  \item{x}{k x 3 matrix containing 3D-coordinates or object of class "mesh3d".

}
  \item{mesh}{triangular surface mesh stored as object of class "mesh3d".
}
  \item{k}{neighbourhood of kd-tree to search - the larger, the slower -
  but the more likely the absolutely closest point is hit.
}
  \item{sign}{logical: if TRUE, signed distances are returned.
}
  \item{cores}{integer: how many cores to use for the search algorithm.
}
  \item{method}{integer: either 0 or 1, if 0 ordinary Euclidean distance
  is used, if 1, the distance suggested by Moshfeghi(1994) is calculated.
}
  \item{\dots}{additional arguments. currently unavailable.
}
}
\details{
  The search for the clostest point is designed as follows:
  Calculate the barycenter of each target face. For each coordinate of
  x, determine the k closest barycenters and calculate the distances to
  the closest point on these faces.
}
\value{
  returns an object of class "mesh3d"
  with:

  \item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
  \item{normals }{4xn matrix containing vertex normals.}
  
  \item{quality }{vector: containing distances to target.}
  \item{it }{4xm matrix containing vertex indices forming triangular
  faces.Only available, when x is a mesh.}

}
\references{ Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating
  Signed Distance Fields From Triangle Meshes. Informatics and
  Mathematical Modelling.

  Moshfeghi M, Ranganath S, Nawyn K. 1994. Three-dimensional elastic matching of volumes IEEE Transactions on Image Processing: A Publication of the IEEE Signal Processing Society 3:128-138.
  
}
\author{
Stefan Schlager}


\seealso{
\code{\link{ply2mesh}}
}
\examples{
data(nose)
out <- closemeshKD(longnose.lm,shortnose.mesh,sign=TRUE)
### show distances - they are very small because
###longnose.lm is scaled to unit centroid size.
hist(out$quality)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
