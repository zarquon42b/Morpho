\name{proc.weight}
\alias{proc.weight}

\title{
calculate weights inverse to the procrustes-distances from the specified observation.
}
\description{
for calculation of a shape model by averaging the observations neighbouring the configuration in question, it is necessary to calculate weights
by similarity.
}
\usage{
proc.weight(data, number, ref,report=TRUE, reg = 0, log = FALSE, mahalanobis = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{array containing landmark configurations

}
  \item{number}{integer: how many of the neighbours are to be involved.
%%     ~~Describe \code{number} here~~
}
  \item{ref}{integer: position in the array that is used as reference.
%%     ~~Describe \code{ref} here~~
}
 \item{report}{logical: require report about name of the reference.
%%     ~~Describe \code{ref} here~~
}
\item{reg}{numeric: regularise mahalanobis distance by adding reg to the
  diagonal of eigenvalues of the covariance matrix.
%%     ~~Describe \code{ref} here~~
}
\item{log}{logical: use the logarithm of the distances.
%%     ~~Describe \code{ref} here~~
}
\item{mahalanobis}{logical: use mahalanobis distance.
%%     ~~Describe \code{ref} here~~
}
}
\details{

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{data }{dataframe containing id, procrustes distance and weight according to the reference}
\item{reference }{returns dimnames(data)[[3]][ref]}
\item{rho.all }{dataframe containing Procrustes distances to references of all observations}
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


## The function is currently defined as
function (data, number, ref) 
{
    col <- 2
    l <- dim(data)[3]
    rho <- 0
    for (i in 1:l) {
        rho[i] <- angle.calc(data[, , ref], data[, , i])$rho
    }
    if (is.null(dimnames(data)[[3]])) {
        id <- as.character(c(1:l))
    }
    else {
        id <- dimnames(data)[[3]]
    }
    data <- data.frame(id, rho)
    dat.sort.i <- data[order(data[, col]), ]
    dat.which <- dat.sort.i[2:(number + 1), ]
    all <- sum(dat.which[, col])
    share <- 0
    for (i in 1:length(dat.which[, col])) {
        share[i] <- dat.which[i, col]/all
    }
    weight <- sort(share, decreasing = T)
    out <- data.frame(dat.which, weight)
    cat(paste("  reference was", id[ref], "\n"))
    return(list(data = out, reference = id[ref], rho.all = data))
  }
}

