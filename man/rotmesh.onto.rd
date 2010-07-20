\name{rotmesh.onto}
\alias{rotmesh.onto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{rotate a mesh

}
\description{    rotates and reflects a mesh onto a set of landmarks
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rotmesh.onto(obj, refmat, tarmat,adnormals = TRUE, scale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ 3D mesh imported with the read.obj fuction
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
%%  ~Describe the value returned
%%  If it is a LIST, use
    \item{obj }{rotated mesh}
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

\seealso{\code{\link{shapelist3d}},\code{\link{warp_obj}} ,\code{\link{trimesh.obj}},\code{\link{read.obj}} ,\code{\link{write.obj}} 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(obj,refmat,tarmat)
{ rot<-rotonto(tarmat,refmat)
  centmesh<-t(apply(obj[which(obj[,1]=="v"),2:4],1,function(x){x - rot$transy}))
  centmeshr<-centmesh\%*\%rot$gamm
  obj[which(obj[,1]=="v"),2:4]<-t(apply(centmeshr,1,function(x){x+rot$trans}))
  if (sign(det(rot$gamm)<0))
  {obj<-conv2backf(obj)}
  return(list(obj=obj,rot=rot))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
