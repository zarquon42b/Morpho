\name{permudist}
\alias{permudist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{performs permutation testing for group differences.
}
\description{
  This function compares the distance between two groupmeans to the
  distances obtained by random assignment of observations to this groups.
}
\usage{
permudist(data, groups, rounds = 1000, which = 1:2, mc.cores = detectCores())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{array or matrix containing data
%%     ~~Describe \code{data} here~~
}
  \item{groups}{factors determining grouping.
%%     ~~Describe \code{groups} here~~
}
  \item{rounds}{number of permutations
%%     ~~Describe \code{rounds} here~~
}
  \item{which}{in case the factor levels are >2 this determins which
    factorlevels to use
%%     ~~Describe \code{which} here~~
}
\item{mc.cores}{integer: determines how many cores to use for the computation. The default
  is autodetect. But in case, it doesn't work as expected cores can be
  set manually.
%%     ~~Describe \code{initproc} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{permudist }{vector containing distances between random group means.}
  \item{dist }{distance between actual group means}
  \item{p.value }{p-value obtained by comparing the actual distance to
    randomly acquired distances.}

  
%% ...
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
