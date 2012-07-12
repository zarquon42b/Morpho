\name{rotmesh.onto}
\alias{rotmesh.onto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{rotate ,scale and translate a mesh based on landmark information.

}
\description{rotates and reflects a mesh onto by calculating the
  transformation from two sets of referenced landmarks.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rotmesh.onto(mesh, refmat, tarmat,adnormals = TRUE, scale = FALSE
, reflection = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{object of class mesh3d.
%%     ~~Describe \code{obj} here~~
}
  \item{refmat}{ k x m matrix with landmarks on the mesh
%%     ~~Describe \code{refmat} here~~
}
  \item{tarmat}{ k x m matrix as target configuration
%%     ~~Describe \code{tarmat} here~~
}
 \item{adnormals}{logical - if TRUE, vertex normals will be updated after rotation.
%%     ~~Describe \code{tarmat} here~~
}
 \item{scale}{ logical: if TRUE the mesh will be scaled according to the size of the target.
%%     ~~Describe \code{tarmat} here~~
}
 \item{reflection}{logical: allow reflection.
    
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
%%  ~Describe the value returned
%%  If it is a LIST, use
    \item{mesh }{rotated mesh}
    \item{yrot }{rotated refmat}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{  Stefan Schlager
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{file2mesh}},\code{\link{warp.mesh}} ,\code{\link{rotonto}},\code{\link{mesh2ply}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
