\name{proc.weight}
\alias{proc.weight}

\title{
calculate weights inverse to the distances from the specified observation.
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
\details{distances of zero will get a weight of 1e12 (this is scaled to
  all weights summing to one), thus weights for observations further
  away are converging to zero.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{data }{dataframe containing id, procrustes/mahalanobis distance and weight according to the reference}
\item{reference }{returns observations' names if available}
\item{rho.all }{dataframe containing distances to references of all observations}
%% ...
}

\examples{

proc <- procSym(gorf.dat)
##get weights for the the four specimen closest to the first observation.
weights <- proc.weight(proc$rotated,4,1)

##estimate the first specimen by weighted neighbour shapes.
estim <- proc$mshape*0;for (i in 1:4){estim <-estim+proc$rotated[,,weights$data$nr[i]]*weights$data$weight[i]}

### visualise
plot(estim,asp=1)## show estimation
points(proc$rotated[,,1],col=3)##show original
}

