\name{mergeMeshes}
\alias{mergeMeshes}
\title{
merge multiple triangular meshes into a single one}
\description{
merge multiple triangular meshes into a single one, preserving color and
vertex normals.
}
\usage{
mergeMeshes(...)
}
\arguments{
  \item{\dots}{triangular meshes of class 'mesh3d' to merge.

  }
}

\value{returns the meshes merged into a single one.

}
\seealso{\code{\link{mesh2ply}, \link{file2mesh}, \link{ply2mesh}}

}
\examples{
data(boneData)
data(nose)
mergedMesh <- mergeMeshes(shortnose.mesh, skull_0144_ch_fe.mesh)
require(rgl)
shade3d(mergedMesh, col=3)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
