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
ignore.stdout = FALSE, sign = FALSE, lmdump=NULL, prodump=NULL)
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
\item{smooth}{logical: rerturn smoothed normals.
%%     ~~Describe \code{readnormals} here~~
}
 \item{ignore.stdout}{logical: supress console output.
%%     ~~Describe \code{readnormals} here~~
}
\item{sign}{logical: request signed distances.
%%     ~~Describe \code{raytol} here~~
}
\item{prodump}{character: name of the dumpfile storing the projected
  points' coordinates (useful when
  using parallel backend).

%%     ~~Describe \code{raytol} here~~
}
\item{lmdump}{character: name of the dumpfile storing data to be projected (useful when
  using parallel backend).

%%     ~~Describe \code{raytol} here~~
}
}

\value{if readnormals = FALSE, a m x 3 matrix containing projected
  points is returned, otherwise a list, where
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
\note{The usage of this function requires the command line tools from
  trimesh-tools (https://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/) installed.
 
}



\seealso{
\code{\link{closemeshKD}}
}
\examples{

data(nose)
\dontrun{repro <- projRead(shortnose.lm,shortnose.mesh)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
