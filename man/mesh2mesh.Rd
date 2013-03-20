\name{mesh2mesh}
\alias{mesh2mesh}
\alias{ray2mesh}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{projects the vertices of a mesh onto the surface of another one.

}
\description{projects the vertices of a mesh onto the surface of another
  one by searching for the closest point (Euclidean distance or along
  vertex normals) on the target by for each vertex.

}
\usage{
mesh2mesh(mesh1, tarmesh, clean = TRUE, cloud = FALSE, sign = FALSE)
ray2mesh(mesh1, tarmesh, tol = 1, angmax = NULL, clean = TRUE, 
    outname = NULL, readback = TRUE, inbound = FALSE, strict = FALSE, 
    ignore.stdout = FALSE, mindist = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh1}{mesh to project. Can be  an object of class "mesh3d" or
    path to an external mesh file (ply, obj, stl).
%%     ~~Describe \code{mesh1} here~~
}
  \item{tarmesh}{mesh to project onto. Can be  an object of class "mesh3d" or
    path to an external mesh file (ply, obj, stl).
}
  \item{clean}{logical: request removing of dumpfiles.
%%     ~~Describe \code{clean} here~~
}
  \item{cloud}{logical: if TRUE and mesh1 is an external mesh file, face
    information will not be read from this file (saves time), as there
    is none present.
%%     ~~Describe \code{clean} here~~
}
 \item{sign}{logical: if TRUE, signed distances are returned.
 }
  \item{angmax}{numeric: maximum angle (in radians) of normals of
    original and hit point are allowed to differ.
  }
  \item{outname}{character: set name of dumpfile used in the process.
  }
  \item{readback}{logical: whether the data is to be read into workspace.
  }
  \item{inbound}{inverse search direction along rays.
  }
   \item{strict}{logical: writes the value 1e12 into vertex quality if no face is hit by the ray.
   }
    \item{ignore.stdout}{suppress command line output.
    }
    \item{mindist}{search in both directions of the ray and use closest points.
  }
}
\details{needs trimesh_project and rayproject from trimesh-tools to be installed \url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/})

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mesh}{projected mesh}

}
\references{ Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating
  Signed Distance Fields From Triangle Meshes. Informatics and
  Mathematical Modelling.

}
\author{Stefan Schlager
}

\seealso{
 \code{\link{ply2mesh}},  \code{\link{closemeshKD}}
}
\examples{
}

