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
 1e-10, mc.cores=detectCores())
}
\arguments{
  \item{data}{array or matrix containing data.
}
  \item{groups}{factors of firs two grouping variables.
}
  \item{subgroups}{factors of the subgrouping.
}
  \item{rounds}{number of requested permutation rounds
}
  \item{scale}{if TRUE: data will be scaled by pooled within group
 covarivance matrix. Otherwise Euclidean distance will be used for
 calculating distances.
}
  \item{tol}{threshold for inverting covariance matrix.
}
\item{mc.cores}{integer: determines how many cores to use for the computation. The default
  is autodetect. But in case, it doesn't work as expected cores can be
  set manually.Parallel processing is disabled on
    Windows due to occasional errors.
}
}
\details{
This function calculates  means of all four subgroups and compares the
 residual vectors of the major grouping variables by angle and distance. 
}
\value{
  \item{angle }{angle between the vectors of the subgroups means}
  \item{dist }{ditances between subgroups }
  \item{meanvec }{matrix containing the means of all four subgroups.}
  \item{permutangles }{vector containing angles from random permutation. }
  \item{permudists }{vector containing distances from random permutation.}
  \item{p.angle }{p-value of angle between residual vectors.}
  \item{p.dist }{p-value of length difference between residual vectors.}
  \item{subdist }{length of residual vectors connecting the subgroups means.}

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
