\name{permudist}
\alias{permudist}
\title{performs permutation testing for group differences.
}
\description{
  This function compares the distance between two groupmeans to the
  distances obtained by random assignment of observations to this groups.
}
\usage{
permudist(data, groups, rounds = 1000, which = 1:2, mc.cores = detectCores())
}
\arguments{
  \item{data}{array or matrix containing data
}
  \item{groups}{factors determining grouping.
}
  \item{rounds}{number of permutations
}
  \item{which}{in case the factor levels are >2 this determins which
    factorlevels to use
}
\item{mc.cores}{integer: determines how many cores to use for the computation. The default
  is autodetect. But in case, it doesn't work as expected cores can be
  set manually. Parallel processing is disabled on
    Windows due to occasional errors.
}
}

\value{
  \item{permudist }{vector containing distances between random group means.}
  \item{dist }{distance between actual group means}
  \item{p.value }{p-value obtained by comparing the actual distance to
    randomly acquired distances.}

  
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
