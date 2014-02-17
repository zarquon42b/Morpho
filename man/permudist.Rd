\name{permudist}
\alias{permudist}
\title{performs permutation testing for group differences.}
\usage{
permudist(data, groups, rounds = 1000, which = NULL)
}
\arguments{
  \item{data}{array or matrix containing data}

  \item{groups}{factors determining grouping.}

  \item{rounds}{number of permutations}

  \item{which}{integer (optional): in case the factor
  levels are > 2 this determins which factorlevels to use}
}
\value{
\item{dist }{distance matrix with distances between actual
group means} \item{p.value }{distance matrix containing
pairwise p-values obtained by comparing the actual distance
to randomly acquired distances}
}
\description{
This function compares the distance between two groupmeans
to the distances obtained by random assignment of
observations to this groups.
}
\examples{
data(boneData)
proc <- procSym(boneLM)
groups <- name2factor(boneLM,which=3)
perm <- permudist(proc$PCscores[,1:10], groups=groups, rounds=10000)

## now we concentrate only on sex dimorphism between Europeans
groups <- name2factor(boneLM,which=3:4)
levels(groups)
perm1 <- permudist(proc$PCscores, groups=groups,which=3:4, rounds=10000)
}

