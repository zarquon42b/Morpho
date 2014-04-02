\name{cutSpace}
\alias{cutSpace}
\title{separate a 3D-pointcloud by a hyperplane
separate a 3D-pointcloud by a hyperplane}
\usage{
cutSpace(pointcloud, v1, v2, v3, upper = TRUE)
}
\arguments{
  \item{pointcloud}{numeric n x 3 matrix}

  \item{v1}{numeric vector of length=3 specifying a point
  on the separating plane}

  \item{v2}{numeric vector of length=3 specifying a point
  on the separating plane}

  \item{v3}{numeric vector of length=3 specifying a point
  on the separating plane}

  \item{upper}{logical specify whether the points above or
  below the plane are to be reported as TRUE.}
}
\value{
logical vector of length n. Reporting for each point if it
is above or below the hyperplane
}
\description{
separate a 3D-pointcloud by a hyperplane separate a
3D-pointcloud by a hyperplane
}
\details{
As above and below are specified by the normal calculated
from \eqn{(v2-v1) \times (v3-v1)}, where \eqn{\times} is
the vector crossproduct. This means the normal points
"upward" when viewed from the positon where v1, v2 and v3
are arranged counter-clockwise. Thus, which side is "up"
depends on the ordering of v1, v2 and v3.
}
\examples{
data(nose)
v1 <- shortnose.lm[1,]
v2 <- shortnose.lm[2,]
v2 <- shortnose.lm[3,]
pointcloud <- vert2points(shortnose.mesh)
upper <- cutSpace(pointcloud, v1, v2, v3)
\dontrun{
points3d(pointcloud[upper,])
}
}

