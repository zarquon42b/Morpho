\name{mesh2ply}
\alias{mesh2ply}
\alias{mesh2obj}
\title{
export mesh objects to disk.
}
\description{
export an object of class mesh3d to a common mesh file.
}
\usage{
mesh2ply(x, filename = "default",col = NULL, writeNormals=FALSE)
mesh2obj(x, filename = "default")
}
\arguments{
  \item{x}{object of class mesh3d - see rgl documentation for further
    details. (At the moment only triangular meshes are supported)
}
  \item{filename}{character: Path/name of the requested output - extension will be added manually.
}
 \item{col}{ Writes color information to ply file. Can be either a single color value or a vector containing a color value for each vertex of the mesh.
}
 \item{writeNormals}{logical: if TRUE, existing normals of a mesh are
   written to file - can slow things down for very large meshes.
}
}

\author{Stefan Schlager
}


\seealso{
\code{\link{ply2mesh}}
}
\examples{
vb<-c(-1.8,-1.8,-1.8,1.0,1.8,-1.8,-1.8,1.0,-1.8,1.8,-1.8,1.0,1.8,
1.8,-1.8,1.0,-1.8,-1.8,1.8,1.0,1.8,
-1.8,1.8,1.0,-1.8,1.8,1.8,1.0,1.8,1.8,1.8,1.0)
it<-c(2,1,3,3,4,2,3,1,5,5,7,3,5,1,2,2,6,5,6,8,7,7,5,6,7,8,4,4,3,7,4,8,6,6,2,4)
vb<-matrix(vb,4,8) ##create vertex matrix
it<-matrix(it,3,12) ## create face matrix
cube<-list(vb=vb,it=it)
class(cube)<-"mesh3d"
shade3d(cube,col=3) ## view the green cube
mesh2ply(cube,filename="cube") # write cube to a file called cube.ply
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

