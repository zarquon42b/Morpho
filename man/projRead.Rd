\name{projRead}
\alias{projRead}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Project points onto the closest point on a mesh

}
\description{project points onto a given surface and return projected points and normals.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
projRead(lm, mesh, readnormals = TRUE, clean = TRUE, smooth = TRUE,
ignore.stdout = FALSE, sign=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lm}{m x 3 matrix containing 3D coordinates.
%%     ~~Describe \code{lm} here~~
}
  \item{mesh}{character: specify path to mesh file.
%%     ~~Describe \code{mesh} here~~
}
  \item{readnormals}{logical: return normals of projected points.
%%     ~~Describe \code{readnormals} here~~
}
  \item{clean}{logical: remove dumpfiles.
%%     ~~Describe \code{readnormals} here~~
}
  \item{smooth}{logical: rerurn smoothed normals.
%%     ~~Describe \code{readnormals} here~~
}
 \item{ignore.stdout}{logical: supress console output.
%%     ~~Describe \code{readnormals} here~~
}
\item{sign}{logical: request signed distances.
%%     ~~Describe \code{raytol} here~~
}
}
\details{
projRead needs the command line tool trimesh_project (http://sourceforge.net/projects/morpho-rpackage/files/) installed.
}
\value{either a m x3 matrix containing projected points or a list, where
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{vb }{3 x m matrix containing projected points}
\item{normals }{3 x m matrix containing normals}
%% ...
}
\references{ Detection of inside/outside uses the algorithm proposed in:

  Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance Fields From Triangle Meshes. Informatics and Mathematical Modelling.
}

\author{
Stefan Schlager}
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
