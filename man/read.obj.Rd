\name{read.obj}
\alias{read.obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Importing a 3D-Mesh
}
\description{
imports a .obj file into R.
}
\usage{
read.obj(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{a obj-file containing information about vertices and facets.
%%     ~~Describe \code{file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{out }{dataframe containing information about vertices and facets}
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

\section{Warning }{be sure that the .obj-file only contains information about vertices and facets. The easy way is to convert your mesh to an obj-file via the free Programm meshlab (http://meshlab.sourceforge.net). }

\seealso{\code{\link{warp_obj}}  ,\code{\link{write.obj}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file) 
{
    out <- read.table(file, sep = " ")
    vn <- (grep("vn", out[, 1]))
    if (length(vn) != 0) {
        out <- out[-vn, ]
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
