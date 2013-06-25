\name{fixLMtps}
\alias{fixLMtps}
\title{estimate missing landmarks
}
\description{Missing landmarks are estimated by deforming a sample average or a weighted estimate of the configurations most similar onto the deficient configuration. The deformation is performed by a Thin-plate-spline interpolation calculated by the available landmarks.
}
\usage{
fixLMtps(data,comp=3,weight=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{array containing landmark data
%%     ~~Describe \code{data} here~~
}
  \item{comp}{integer: select how many of the closest observations are to be taken to calculate an initial estimate.
%%     ~~Describe \code{data} here~~
}
  \item{weight}{logical: requests the calculation of an estimate based
    on the procrustes distance. Otherwise the sample's consensus is used
    as reference.
%%     ~~Describe \code{data} here~~
}
}
\details{This function tries to estimate missing landmark data by mapping weighted averages from complete datasets onto the missing specimen. The weights are the inverted Procrustes (see  \code{\link{proc.weight}}) distances between the 'comp' closest specimen (using the available landmark configuration).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{out }{array containing all data, including fixed configurations  - same order as input}
\item{mshape }{meanshape  - calculated from complete datasets.}
\item{checklist }{list containing information about missing landmarks}
\item{check }{vector containing position of observations in data where at least one missing coordinate was found }
%% ...
}
\references{Bookstein FL. 1989. Principal Warps: Thin-plate splines and the decomposition of deformations IEEE Transactions on pattern analysis and machine intelligence 11.

}
\author{Stefan Schlager
%%  ~~who you are~~
}
\note{Be aware that these estimates might be grossly wrong when the
  missing landmark is quite far off the rest of the landmarks (due to
  the radial basis function used in the Thin-plate spline interpolation.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{proc.weight}},  \code{\link{tps3d}}
}
\examples{
library(shapes)
data <- gorf.dat
### set first landmark of first specimen to NA
data[1,,1] <- NA
repair <- fixLMtps(data,comp=5)
### view difference between estimated and actual landmark
plot(repair$out[,,1],asp=1,pch=21,cex=0.7,col=2)#estimated landmark
points(gorf.dat[,,1],col=3,pch=20)#actual landmark



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
