\name{unrefVertex}
\alias{unrefVertex}
\alias{rmVertex}
\alias{vert2points}
\alias{rmUnrefVertex}
\title{some little helpers for vertex operations on triangular meshes
}
\description{
remove (unreferenced vertices from triangular meshes)}
\usage{
unrefVertex(mesh)
rmVertex(mesh,index)
vert2points(mesh)
rmUnrefVertex(mesh)

}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d".
  }
  \item{index}{vector containing indices of vertices to be removed.
}
}
\value{unrefVertex: vector with indices of unreferenced vertices.
  
  rmVertex: returns mesh with specified vertices removed and faces and
  normals updated.
  
  vert2points: k x 3 matrix containing vertex coordinates.
  
  rmUnrefVertex: mesh with unreferenced vertices removed.
  
}

\author{
Stefan Schlager}



\seealso{
 \code{\link{ply2mesh}}, \code{\link{file2mesh}} 
}
\examples{
data(nose)
testmesh <- rmVertex(shortnose.mesh,1:50) ## remove first 50 vertices
shade3d(testmesh,col=3) ### view result
testmesh$vb <- cbind(testmesh$vb,shortnose.mesh$vb[,1:50]) ## add some unreferenced vertices
points3d(vert2points(testmesh),col=2)## see the vertices in the holes?

cleanmesh <- rmUnrefVertex(testmesh)## remove those lonely vertices!
rgl.pop()
points3d(vert2points(cleanmesh),col=2) ### now the holes are empty!!


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
