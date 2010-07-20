\name{write.obj}
\alias{write.obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{exporting a 3D mesh into an obj-file
}
\description{
writes an imported (and manipulated) .obj mesh into a file
}
\usage{
write.obj(obj, filename = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{via read.obj imported 3D mesh.
%%     ~~Describe \code{obj} here~~
}
  \item{filename}{character string for the filename of the generated .obj-file.
%%     ~~Describe \code{filname} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{warp_obj}} ,\code{\link{read.obj}} }
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(obj,filname="default")
{write.table(obj,file=paste(filename,".obj"),quote=F,row.names = FALSE,col.names = FALSE)}
}
