\name{regdist.raw}
\alias{regdist.raw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{correlation between shape space and tangent space

}
\description{performs a full procrustes superimposition of landmark data and calculates the correaltion between tangent and shape space.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regdist.raw(dataarray, plot = TRUE, main = "", rho = "angle", dist.mat.out = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataarray}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
%%     ~~Describe \code{dataarray} here~~
}
  \item{plot}{Logical: whether to plot the distances between observations.
%%     ~~Describe \code{plot} here~~
}
  \item{main}{character string: Title of the plot.
%%     ~~Describe \code{main} here~~
}
  \item{rho}{chose how to calculate distances in shape space. Options: "riemdist"=Riemannian distance (function from the shapes package-takes along time to calculate), "angle"=calculates the angle between shape vectors, "sindist"=sinus of length of residual vector between shape vectors.
%%     ~~Describe \code{rho} here~~
}
  \item{dist.mat.out}{Logical: If TRUE, output will contain distance matrices.
%%     ~~Describe \code{dist.mat.out} here~~
}
%%  \item{reflect}{Logical: indicating if reflection is required for Procrustes superimposition.
%%     ~~Describe \code{reflect} here~~
%%}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{cor}{correlation coefficient between distances in shape space and tangent space}
\item{procSS }{Procrustes Sums of Squares (of full procrustes distance.}
\item{tanSS }{Tangent Sums of Squares}
\item{rhoSS }{Procrustes Sums of Squares (of angle) }
\item{euc.dist }{distance matrix of euclidean distance in Tangent space}
\item{proc.dist }{distance matrix of Procrustes distance in Shape space}
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
Stefan Schlager
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{regdist}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dataarray, plot = TRUE, main = "", rho = "angle", dist.mat.out = FALSE, 
    reflect = FALSE) 
{
    proc <- procSym(dataarray, reflect = reflect)
    x <- proc$rotated
    n <- dim(x)[3]
    m <- dim(x)[2]
    k <- dim(x)[1]
    y <- proc$orpdata
    qm <- dist(t(matrix(x, k * m, n)))
    procdis <- sum(qm^2)/n
    procdistmat <- matrix(NA, n, n)
    for (i in 1:n) {
        for (j in 1:n) if (rho == "riemdist") {
            {
                procdistmat[i, j] <- riemdist(x[, , i], x[, , 
                  j])
            }
        }
        else if (rho == "angle") {
            {
                procdistmat[i, j] <- angle.calc(x[, , i], x[, 
                  , j])$rho
            }
        }
    }
    if (rho == "sindist") {
        procvec <- asin(qm)
    }
    else {
        procvec <- as.dist(procdistmat)
    }
    procdis2 <- sum(procvec^2)/n
    em <- dist(t(matrix(y, k * m, n)))
    euvec <- (em)
    eudis <- sum(euvec^2)/n
    correlation <- cor(euvec, procvec)^2
    if (plot == TRUE) {
        plot(euvec, procvec, asp = 1, xlab = "euclid. dist. in tangentspace", 
            ylab = paste("rho as", rho), main = main)
        abline(0, 1, col = "grey50")
    }
    if (dist.mat.out == TRUE) {
        return(list(cor = correlation, procSS = procdis, tanSS = eudis, 
            rhoSS = procdis2, euc.dist = em, proc.dist = procvec))
    }
    else {
        return(list(cor = correlation, procSS = procdis, tanSS = eudis, 
            rhoSS = procdis2))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
