\name{read.lmdta}
\alias{read.lmdta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{read dta files
}
\description{
reads .dta files created by the software Landmark http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\usage{
read.lmdta(file = "x")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{a dta file 
%%     ~~Describe \code{file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{arr}{array containing landmarks dimnames will be Information of ID and landmark names specified in Landmark}
 \item{info }{Information extracted from the header of the dta file}
 \item{idnames }{character vector containing the names of the individuals as specified in the dta file}
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
    x <- file
    A <- readLines(x)
    em <- which(A == "")
    idnames <- A[c((em[1] + 1):(em[2] - 1))]
    info <- strsplit(A[3], " ")[[1]]
    n2 <- nchar(info[2]) - 1
    nspeci <- as.numeric(substr(info[2], 1L, n2))
    ndim <- as.numeric(substr(info[6], 5, nchar(info[6])))
    nlms <- as.numeric(info[3])/ndim
    eot <- em[2]
    B <- as.matrix(read.table(x, skip = eot), na.strings = as.numeric(info[5]))
    tt <- array(t(B), dim = c(ndim, nlms, nspeci))
    arr <- array(NA, dim = c(nlms, ndim, nspeci))
    for (i in 1:nspeci) {
        arr[, , i] <- t(tt[, , i])
    }
    dimnames(arr)[[3]] <- idnames
    return(list(arr = arr, info = info, idnames = idnames))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
