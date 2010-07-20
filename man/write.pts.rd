\name{write.pts}
\alias{write.pts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
exports a matrix containing landmarks into .pts format
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.pts(X, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{k x m matrix containing landmark configuration
%%     ~~Describe \code{X} here~~
}
  \item{filename}{ filename - extension .pts will be added automatically
%%     ~~Describe \code{filename} here~~
}
}
\details{you can import the information into the program landmarks available at http://graphics.idav.ucdavis.edu/research/EvoMorph
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
\author{  Stefan Schlager
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
function (X, filename) 
{
    filename <- paste(filename, ".pts", sep = "")
    k <- dim(X)[1]
    m <- dim(X)[2]
    X <- as.matrix(X)
    a0 <- paste("S", 0, c(0:(k - 1)), sep = "")
    all.frame <- cbind(a0, X)
    cat("Version 1.0\n", file = filename)
    cat(paste(k, "\n", sep = ""), file = filename, append = T)
    write(t(all.frame), file = filename, sep = " ", ncolumns = m + 
        1, append = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
