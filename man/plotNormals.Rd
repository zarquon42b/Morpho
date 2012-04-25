\name{plotNormals}
\alias{plotNormals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plots the normals of a triangular surface mesh.
}
\description{
  visualises the vertex normals of a triangular surface mesh of class
  mesh3d. If no normals are contained, they are computed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotNormals(x, long = 1, lwd = 1, col = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class "mesh3d"
}
  \item{long}{length of the normals (default is 1)
%%     ~~Describe \code{long} here~~
}
  \item{lwd}{width of the normals
%%     ~~Describe \code{lwd} here~~
}
  \item{col}{color of the normals
%%     ~~Describe \code{col} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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

data(nose)
shade3d(shortnose.mesh,col=3)
plotNormals(shortnose.mesh,col=4,long=0.01)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
