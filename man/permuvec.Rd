\name{permuvec}
\alias{permuvec}
\title{ perfom permutation testing on angles and distances between
  subgroups of two major groups.
}
\description{

 perform permutation test on length and angle of the vectors connecting
 the subgroup means of two groups: e.g. compare if length and angle
 between sex related differences in two populations differ significantly.
}
\usage{
permuvec(data, groups, subgroups, rounds = 10000, scale = TRUE, tol =
 1e-10, mc.cores=detectCores()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{array or matrix containing data.
%%     ~~Describe \code{data} here~~
}
  \item{groups}{factors of firs two grouping variables.
%%     ~~Describe \code{groups} here~~
}
  \item{subgroups}{factors of the subgrouping.
%%     ~~Describe \code{subgroups} here~~
}
  \item{rounds}{number of requested permutation rounds
%%     ~~Describe \code{rounds} here~~
}
  \item{scale}{if TRUE: data will be scaled by pooled within group
 covarivance matrix. Otherwise Euclidean distance will be used for
 calculating distances.
%%     ~~Describe \code{scale} here~~
}
  \item{tol}{threshold for inverting covariance matrix.
%%     ~~Describe \code{tol} here~~
}
\item{mc.cores}{integer: determines how many cores to use for the computation. The default
  is autodetect. But in case, it doesn't work as expected cores can be
  set manually.
%%     ~~Describe \code{initproc} here~~
}
}
\details{
This function calculates  means of all four subgroups and compares the
 residual vectors of the major grouping variables by angle and distance. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{angle }{angle between the vectors of the subgroups means}
  \item{dist }{ditances between subgroups }
  \item{meanvec }{matrix containing the means of all four subgroups.}
  \item{permutangles }{vector containing angles from random permutation. }
  \item{permudists }{vector containing distances from random permutation.}
  \item{p.angle }{p-value of angle between residual vectors.}
  \item{p.dist }{p-value of length difference between residual vectors.}
  \item{subdist }{length of residual vectors connecting the subgroups means.}

%% ...
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
