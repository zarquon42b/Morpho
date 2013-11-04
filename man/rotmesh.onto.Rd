\name{rotmesh.onto}
\alias{rotmesh.onto}
\title{rotate ,scale and translate a mesh based on landmark information.

}
\description{rotates and reflects a mesh onto by calculating the
  transformation from two sets of referenced landmarks.
}
\usage{
rotmesh.onto(mesh, refmat, tarmat, adnormals = FALSE, scale = FALSE
, reflection = FALSE)
}
\arguments{
  \item{mesh}{object of class mesh3d.
}
  \item{refmat}{ k x m matrix with landmarks on the mesh
}
  \item{tarmat}{ k x m matrix as target configuration
}
 \item{adnormals}{logical - if TRUE, vertex normals will be recomputed
   after rotation. If \code{mesh} has normals and adnormals=FALSE, the existing
   normals are rotated by the same rotation matrix as the mesh's vertices.
}
 \item{scale}{ logical: if TRUE the mesh will be scaled according to the size of the target.
}
 \item{reflection}{logical: allow reflection.
    
  }
}

\value{ 
    \item{mesh }{rotated mesh}
    \item{yrot }{rotated refmat}
}
\author{  Stefan Schlager
}



\seealso{\code{\link{file2mesh}},\code{\link{warp.mesh}} ,\code{\link{rotonto}},\code{\link{mesh2ply}}
}
\examples{
require(rgl)
data(boneData)
## rotate, translate and scale the mesh belonging to the first specimen
## onto the landmark configuration of the 10th specimen
rotmesh <- rotmesh.onto(skull_0144_ch_fe.mesh,boneLM[,,1],
                        boneLM[,,10], scale=TRUE)
\dontrun{
## render rotated mesh and landmarks
shade3d(rotmesh$mesh, col=2, specular=1)
spheres3d(boneLM[,,1])
## render original mesh
shade3d(skull_0144_ch_fe.mesh, col=3, specular=1)
spheres3d(boneLM[,,10])
}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
