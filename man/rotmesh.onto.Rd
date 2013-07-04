\name{rotmesh.onto}
\alias{rotmesh.onto}
\title{rotate ,scale and translate a mesh based on landmark information.

}
\description{rotates and reflects a mesh onto by calculating the
  transformation from two sets of referenced landmarks.
}
\usage{
rotmesh.onto(mesh, refmat, tarmat,adnormals = TRUE, scale = FALSE
, reflection = FALSE)
}
\arguments{
  \item{mesh}{object of class mesh3d.
}
  \item{refmat}{ k x m matrix with landmarks on the mesh
}
  \item{tarmat}{ k x m matrix as target configuration
}
 \item{adnormals}{logical - if TRUE, vertex normals will be updated after rotation.
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
