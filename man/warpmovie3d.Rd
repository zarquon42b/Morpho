\name{warpmovie3d}
\alias{warpmovie3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a sequence of images showing predefined steps of warping two 
meshes into each other}
\description{given two meshes with the same amount of vertices and faces
  (e.g a mesh and its warped counterpart), the starting mesh will be
  subsequently transformed into the final mesh by splitting the
  differences into a predefined set of steps.
}
\usage{
warpmovie3d(x, y, n, col = skin1, palindrome = FALSE, folder = NULL, movie = "warpmovie", add = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{mesh to start with (object of class mesh3d)
%%     ~~Describe \code{x} here~~
}
  \item{y}{resulting mesh (object of class mesh3d), having the same
    amount of vertices and faces than the starting mesh
%%     ~~Describe \code{y} here~~
}
  \item{n}{integer: amount of intermediate steps.
%%     ~~Describe \code{n} here~~
}
  \item{col}{color of the mesh
%%     ~~Describe \code{col} here~~
}
  \item{palindrome}{logical: if TRUE, the procedure will go forth and back.
%%     ~~Describe \code{palindrome} here~~
}
  \item{folder}{character: output folder for created images (optional)
%%     ~~Describe \code{folder} here~~
}
  \item{movie}{character: name of teh output files
%%     ~~Describe \code{movie} here~~
}
  \item{add}{logical: if TRUE, the movie will be added to the focussed rgl-windows.
%%     ~~Describe \code{add} here~~
}
}
\details{
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
    \code{\link{ply2mesh},\link{file2mesh},\link{mesh2ply},\link{warp.mesh}}
  }

  \examples{
data(nose)##load data
##warp a mesh onto another landmark configuration:
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
warpmovie3d(shortnose.mesh,warpnose.long,n=15)## create 15 images.
### the images are now stored in your current working directory and can
### be concatenated to a gif using an external program such as imagemagick.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
