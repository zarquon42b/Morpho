\name{plotAtlas}
\alias{plotAtlas}

\title{
visualize an atlas defined by createAtlas}
\description{
visualize an atlas defined by createAtlas
}
\usage{
plotAtlas(atlas, radius = 1, meshcol = "white", add = TRUE)
}
\arguments{
  \item{atlas}{object of class atlas created by \code{\link{createAtlas}}.
}
  \item{radius}{radius of the spheres representing the fix landmarks.
}
  \item{meshcol}{color to render the atlas mesh
}
  \item{add}{logical: if TRUE, a new rgl window is opened.
}
}
\details{additonally a plot with a legend will open where coloring of
  the 3D-spheres is specified.
}

\seealso{\code{\link{placePatch}, \link{createAtlas}}
}
\examples{
data(nose)
atlas <- createAtlas(shortnose.mesh, landmarks =
           shortnose.lm[c(1:5,20:21),], patch=shortnose.lm[-c(1:5,20:21),])
plotAtlas(atlas)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
