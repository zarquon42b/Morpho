\name{adnormals}
\alias{adnormals}
\title{Compute normals of a triangular mesh
}
\description{
Compute vertex normals of a triangular mesh of class "mesh3d"}
\usage{
adnormals(x, angle = TRUE)
}
\arguments{
  \item{x}{triangular mesh of class "mesh3d"
%%     ~~Describe \code{x} here~~
}
  \item{angle}{logical: if TRUE, angle weighted normals are used.
%%     ~~Describe \code{angle} here~~
}
}

\value{returns mesh with updated vertex normals

}
\references{Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance Fields From Triangle Meshes. Informatics and Mathematical Modelling, . 

}
\author{
Stefan Schlager
}
\note{
only supports triangular meshes}


\seealso{
 \code{\link{ply2mesh}}
}
\examples{
data(nose)
shortnose.mesh$normals <- NULL ##remove normals
shade3d(shortnose.mesh,col=3)##render
shortnose.mesh <- adnormals(shortnose.mesh)
rgl.clear()
shade3d(shortnose.mesh,col=3)##smoothly rendered now
  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
