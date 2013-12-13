\name{mesh2mesh}
\alias{mesh2mesh}
\alias{ray2mesh}
\title{projects the vertices of a mesh onto the surface of another one.}
\usage{
mesh2mesh(mesh1, tarmesh, clean = TRUE, cloud = FALSE, sign = FALSE)

ray2mesh(mesh1, tarmesh, tol = 1, angmax = NULL, clean = TRUE,
  outname = NULL, readback = TRUE, inbound = FALSE, strict = FALSE,
  ignore.stdout = FALSE, mindist = FALSE, targetdump = "target",
  refdump = "reference")
}
\arguments{
  \item{mesh1}{mesh to project. Can be an object of class
  "mesh3d" or path to an external mesh file (ply, obj,
  stl).}

  \item{tarmesh}{mesh to project onto. Can be an object of
  class "mesh3d" or path to an external mesh file (ply,
  obj, stl).}

  \item{clean}{logical: request removing of dumpfiles.}

  \item{cloud}{logical: if TRUE and mesh1 is an external
  mesh file, face information will not be read from this
  file (saves time), as there is none present.}

  \item{sign}{logical: if TRUE, signed distances are
  returned.}

  \item{tol}{numeric: maximum distance to search along ray,
  closest Euclidean distance will be used, if tol is
  exceeded.}

  \item{angmax}{numeric: maximum angle (in radians) of
  normals of original and hit point are allowed to differ.}

  \item{outname}{character: set name of dumpfile used in
  the process.}

  \item{readback}{logical: whether the data is to be read
  into workspace.}

  \item{inbound}{inverse search direction along rays.}

  \item{strict}{logical: writes the value 1e12 into vertex
  quality if no face is hit by the ray.}

  \item{ignore.stdout}{suppress command line output.}

  \item{mindist}{search in both directions of the ray and
  use closest points.}

  \item{targetdump}{specify name of dumped target mesh}

  \item{refdump}{specify name of dumped reference mesh}
}
\value{
returns projected mesh
}
\description{
projects the vertices of a mesh onto the surface of another
one.
}
\details{
projects the vertices of a mesh onto the surface of another
one by searching for the closest point (Euclidean distance
or along vertex normals) on the target by for each vertex.

needs trimesh_project and rayproject from trimesh-tools to
be installed (
\url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/})
}
\author{
Stefan Schlager
}
\references{
Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating
Signed Distance Fields From Triangle Meshes. Informatics
and Mathematical Modelling.
}
\seealso{
\code{\link{ply2mesh}}, \code{\link{closemeshKD}}
}

