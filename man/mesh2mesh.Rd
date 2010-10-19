\name{mesh2mesh}
\alias{mesh2mesh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{projects the vertices of a mesh onto the surface of another one.

}
\description{projects the vertices of a mesh onto the surface of another one.

}
\usage{
mesh2mesh(mesh1, tarmesh, clean = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh1}{mesh to project (object of class "mesh3d")
%%     ~~Describe \code{mesh1} here~~
}
  \item{tarmesh}{mesh to project onto (object of class "mesh3d")
%%     ~~Describe \code{tarmesh} here~~
}
  \item{clean}{logical: request removing of dumpfiles.
%%     ~~Describe \code{clean} here~~
}
}
\details{needs trimesh_project to be installed (http://sourceforge.net/projects/morpho-rpackage/files/)

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mesh}{projected mesh}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

