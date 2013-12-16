\name{asymPermute}
\alias{asymPermute}
\title{Assess differences in amount and direction of asymmetric variation}
\usage{
  asymPermute(x, groups, rounds = 1000, which = 1:2,
    mc.cores = detectCores())
}
\arguments{
  \item{x}{object of class symproc result from calling
  \code{\link{procSym}} with \code{pairedLM} specified}

  \item{groups}{factors determining grouping.}

  \item{rounds}{number of permutations}

  \item{which}{in case the factor levels are >2 this
  determins which factorlevels to use}

  \item{mc.cores}{integer: determines how many cores to use
  for the computation. The default is autodetect. But in
  case, it doesn't work as expected cores can be set
  manually. Parallel processing is disabled on Windows due
  to occasional errors.}
}
\value{
  \item{dist }{difference between vector lengths of group
  means} \item{angle }{angle between vectors of group
  specific asymmetric deviation} \item{means }{actual group
  averages} \item{p.dist }{p-value obtained by comparing
  the actual distance to randomly acquired distances}
  \item{p.angle }{p-value obtained by comparing the actual
  angle to randomly acquired angles} \item{permudist
  }{vector containing differences between random group
  means' vector lenghts} \item{permuangle }{vector
  containing angles between random group means' vectors}
}
\description{
  Assess differences in amount and direction of asymmetric
  variation
}

