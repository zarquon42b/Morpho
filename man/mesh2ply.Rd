\name{mesh2ply}
\alias{mesh2ply}
\alias{mesh2obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
export mesh objects to disk.
}
\description{
export an object of class mesh3d to a common mesh file.
}
\usage{
mesh2ply(x, filename = "default",col = NULL, writeNormals=FALSE)
mesh2obj(x, filename = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class mesh3d - see rgl documentation for further
    details. (At the moment only triangular meshes are supported)
%%     ~~Describe \code{x} here~~
}
  \item{filename}{character: Path/name of the requested output - extension will be added manually.
%%     ~~Describe \code{filename} here~~
}
 \item{col}{ Writes color information to ply file. Can be either a single color value or a vector containing a color value for each vertex of the mesh.
%%     ~~Describe \code{x} here~~
}
 \item{writeNormals}{logical: if TRUE, existing normals of a mesh are
   written to file - can slow things down for very large meshes.
%%     ~~Describe \code{x} here~~
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
\author{Stefan Schlager
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
vb<-c(-1.8,-1.8,-1.8,1.0,1.8,-1.8,-1.8,1.0,-1.8,1.8,-1.8,1.0,1.8,1.8,-1.8,1.0,-1.8,-1.8,1.8,1.0,1.8,-1.8,1.8,1.0,-1.8,1.8,1.8,1.0,1.8,1.8,1.8,1.0)
it<-c(2,1,3,3,4,2,3,1,5,5,7,3,5,1,2,2,6,5,6,8,7,7,5,6,7,8,4,4,3,7,4,8,6,6,2,4)
vb<-matrix(vb,4,8) ##create vertex matrix
it<-matrix(it,3,12) ## create face matrix
cube<-list(vb=vb,it=it)
class(cube)<-"mesh3d"
shade3d(cube,col=3) ## view the green cube
mesh2ply(cube,filename="cube") # write cube to a file called cube.ply
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

