\name{ply2mesh}
\alias{ply2mesh}
\alias{obj2mesh}
\alias{file2mesh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import 3D surface mesh files
}
\description{
imports 3D mesh files and store them as an R .object of class mesh3d
}
\usage{
ply2mesh(filename, adnormals = TRUE)
obj2mesh(filename, adnormals = TRUE)
file2mesh(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{character: path to file
%%     ~~Describe \code{filename} here~~
}
  \item{adnormals}{Logical: If the file does not contain normal information, they will be calculated in R: Can take some time.
%%     ~~Describe \code{adnormals} here~~
}
}
\details{
file2mesh needs the command line tool ply2ascii installed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{mesh }{list of class mesh3d - see rgl manual for further details}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
