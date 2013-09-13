\name{createAtlas}
\alias{createAtlas}

\title{
Create an atlas needed in placePatch
}
\description{
Create an atlas needed in placePatch
}
\usage{
createAtlas(mesh, landmarks, patch, corrCuves = NULL, patchCurves = NULL)
}
\arguments{
  \item{mesh}{triangular mesh representing the atlas' surface
}
  \item{landmarks}{matrix containing landmarks defined on the atlas, as
    well as on each specimen in the corresponding sample.
}
  \item{patch}{matrix containing semi-landmarks to be projected onto
    each specimen in the corresponding sample.
}
  \item{corrCuves}{integer vector specifiyng the rowindices of
    \code{landmarks} to be curves
    defined on the atlas AND each specimen.
}
  \item{patchCurves}{integer vector specifiyng the rowindices of
    \code{patch} to be curves only
    defined on the atlas.
  }
}

\value{Returns a list of class "atlas".
  Its content is corresponding to argument names.
}

\note{This is a helper function of \code{\link{placePatch}}.
}

\seealso{\code{\link{placePatch}, \link{plotAtlas}}
}
\examples{
data(nose)
atlas <- createAtlas(shortnose.mesh, landmarks =
            shortnose.lm[c(1:5,20:21),], patch=shortnose.lm[-c(1:5,20:21),])
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
