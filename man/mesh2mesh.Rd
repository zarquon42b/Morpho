\name{mesh2mesh}
\alias{mesh2mesh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{projects the vertices of a mesh onto the surface of another one.

}
\description{projects the vertices of a mesh onto the surface of another
  one by searching for the closest point on the target by for each vertex.

}
\usage{
mesh2mesh(mesh1, tarmesh, clean = TRUE, cloud = FALSE)
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
}
\details{needs trimesh_project to be installed (http://sourceforge.net/projects/morpho-rpackage/files/)

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mesh}{projected mesh}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}

