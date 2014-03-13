\name{quad2trimesh}
\alias{quad2trimesh}
\title{converts a mesh containing quadrangular faces into one only consisting of triangles}
\usage{
quad2trimesh(mesh)
}
\arguments{
  \item{mesh}{object of class "mesh3d"}
}
\value{
triangular mesh with updated normals
}
\description{
converts a mesh containing quadrangular faces into one only
consisting of triangles
}
\examples{
Sigma <- diag(3:1) #create a 3D-covariance matrix
require(rgl)
quadmesh <- ellipse3d(Sigma)##create quadmesh
trimesh <- quad2trimesh(quadmesh)# convert to trimesh
}

