\name{asymPermute}
\alias{asymPermute}
\title{Assess differences in amount and direction of asymmetric variation}
\usage{
asymPermute(x, groups, rounds = 1000, which = NULL)
}
\arguments{
  \item{x}{object of class symproc result from calling
  \code{\link{procSym}} with \code{pairedLM} specified}

  \item{groups}{factors determining grouping.}

  \item{rounds}{number of permutations}

  \item{which}{select which factorlevels to use, if NULL,
  all pairwise differences will be assessed after shuffling
  pooled data.}
}
\value{
\item{dist }{difference between vector lengths of group
means} \item{angle }{angle between vectors of group
specific asymmetric deviation} \item{means }{actual group
averages} \item{p.dist }{p-value obtained by comparing the
actual distance to randomly acquired distances}
\item{p.angle }{p-value obtained by comparing the actual
angle to randomly acquired angles} \item{permudist }{vector
containing differences between random group means' vector
lenghts} \item{permuangle }{vector containing angles
between random group means' vectors} \item{groupmeans}{
array with asymmetric displacement per group}
\item{levels}{ character vector containing the factors
used}
}
\description{
Assess differences in amount and direction of asymmetric
variation
}

