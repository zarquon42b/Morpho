\name{projRead}
\alias{projRead}
\title{Project points onto the closest point on a mesh}
\usage{
  projRead(lm, mesh, readnormals = TRUE, clean = TRUE,
    smooth = TRUE, ignore.stdout = FALSE, sign = FALSE,
    lmdump = NULL, prodump = NULL)
}
\arguments{
  \item{lm}{m x 3 matrix containing 3D coordinates.}

  \item{mesh}{character: specify path to mesh file.}

  \item{readnormals}{logical: return normals of projected
  points.}

  \item{clean}{logical: remove dumpfiles.}

  \item{smooth}{logical: rerturn smoothed normals.}

  \item{ignore.stdout}{logical: supress console messages
  from system calls.}

  \item{sign}{logical: request signed distances.}

  \item{prodump}{character: name of the dumpfile storing
  the projected points' coordinates (useful when using
  parallel backend).}

  \item{lmdump}{character: name of the dumpfile storing
  data to be projected (useful when using parallel
  backend).}
}
\value{
  if readnormals = FALSE, a m x 3 matrix containing
  projected points is returned, otherwise a list, where
  \item{vb }{3 x m matrix containing projected points}
  \item{normals }{3 x m matrix containing normals}
}
\description{
  project points onto a given surface and return projected
  points and normals.
}
\note{
  The usage of this function requires the command line
  tools from trimesh-tools
  (https://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/)
  installed.
}
\examples{
data(nose)
\dontrun{
repro <- projRead(shortnose.lm,shortnose.mesh)
}
}
\author{
  Stefan Schlager
}
\references{
  Detection of inside/outside uses the algorithm proposed
  in:

  Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating
  Signed Distance Fields From Triangle Meshes. Informatics
  and Mathematical Modelling.
}
\seealso{
  \code{\link{closemeshKD}}
}
\keyword{~kwd1}
\keyword{~kwd2}

