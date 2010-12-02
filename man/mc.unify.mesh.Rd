\name{mc.unify.mesh}
\alias{mc.unify.mesh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{creates meshes with corresponding vertices.
}
\description{Warps one mesh into the vicinity of another one and projects the vertices of the reference one on the target.

}
\usage{
mc.unify.mesh(mesh1, mesh2, data1, data2, ray = FALSE, tol = NULL, updateNormals = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh1}{object of class "mesh3d": reference mesh
%%     ~~Describe \code{mesh1} here~~
}
  \item{mesh2}{object of class "mesh3d": target mesh
%%     ~~Describe \code{mesh2} here~~
}
\item{ray}{logical: requests search of next vertex along normals
%%     ~~Describe \code{mesh2} here~~
}
\item{tol}{numeric: threshold for maximal distance along normals if ray=TRUE (in both directions of the normal), if nothing is hit within this threshold, the closest point on the target will be used.
%%     ~~Describe \code{mesh2} here~~
}
  \item{data1}{k x 3 matrix conataining reference points (landmarks) on mesh1.
%%     ~~Describe \code{data1} here~~
}
  \item{data2}{k x 3 matrix conataining reference points (landmarks) on mesh2.
%%     ~~Describe \code{data2} here~~
}
  \item{updateNormals}{logical: requests vertex normals to be updated after unification.
%%     ~~Describe \code{updateNormals} here~~
}
}
\details{needs comand line programs "trimesh_project" and "rayproject" to be installed (http://sourceforge.net/projects/morpho-rpackage/files/)

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{unimesh }{a (copy) of mesh2 with vertices corresponding to those of mesh1}
 \item{warp.mesh }{mesh1 warped onto mesh2}

}
\references{

}
\author{Stefan Schlager

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
