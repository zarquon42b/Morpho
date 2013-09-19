\name{plotAtlas}
\alias{plotAtlas}

\title{
visualize an atlas defined by createAtlas}
\description{
visualize an atlas defined by createAtlas
}
\usage{
plotAtlas(atlas, radius=1, alpha=1, render=c("w", "s"), point=c("s", "p"),
          meshcol="white", add=TRUE, legend=TRUE)
}
\arguments{
  \item{atlas}{object of class atlas created by \code{\link{createAtlas}}.
}
  \item{radius}{radius of the spheres representing the fix landmarks.
  }
  \item{alpha}{value between 0 and 1. Sets transparency of mesh 1=opaque
    0= fully transparent.
  }
  \item{render}{if \code{render="w"}, a wireframe will be drawn, if \code{render="s"}, the
    mesh will be shaded.
  }
  \item{point}{how to render landmarks. "s"=spheres, "p"=points.
  }
  \item{meshcol}{color to render the atlas mesh
  }
  \item{add}{logical: if TRUE, a new rgl window is opened.
  }
  \item{legend }{logical: request plot of legend specifying landmark
    coloring.
  }
}
\details{If \code{legend=TRUE}, a plot with a legend will open where coloring of
  the 3D-spheres is specified.
}
\value{
  returns invisible vector containing \code{rgl.id} of rendered objects.
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
