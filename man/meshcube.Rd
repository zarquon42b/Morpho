\name{meshcube}
\alias{meshcube}
\title{
creates a bounding box for a mesh }
\description{creates a bounding box for a mesh
}
\usage{
meshcube(x)
}
\arguments{
  \item{x}{object of class 'mesh3d'
}
}

\value{returns a 8 x 3 matrix with the coordinates of the corners of the bounding box.
}

\examples{
data(boneData)
mc <- meshcube(skull_0144_ch_fe.mesh)
spheres3d(mc)
wire3d(skull_0144_ch_fe.mesh)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
