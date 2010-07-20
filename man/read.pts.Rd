\name{read.pts}
\alias{read.pts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reads pts files
}
\description{
reads Landmark data exported from the software Landmark from http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\usage{
read.pts(file = "x")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{pts file 
%%     ~~Describe \code{file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{matrix}{A Matrix containing landmark information rownames will be the names given to the landmarks in Landmark}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file = "x") 
{
    pts <- read.table(file, skip = 2)
    rownames(pts) <- pts[, 1]
    pts <- pts[, 2:4]
    return(pts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
