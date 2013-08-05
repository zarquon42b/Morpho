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
\examples{
data(boneData)
proc <- procSym(boneLM)
groups <- name2factor(boneLM,which=3)
perm <- permudist(proc$PCscores, groups=groups, mc.cores=2)
#visualize results
hist(perm$permudist, xlim=c(0,0.1),main="measured vs. random distances",
     xlab="distances")#random distances
points(perm$dist,50,col=2,pch=19)#actual distance
text(perm$dist,70,label=paste("actual distance\n(p=",perm$p.value,")"))

## now we concentrate only on sex dimorphism between Europeans
groups <- name2factor(boneLM,which=3:4)
levels(groups)
perm1 <- permudist(proc$PCscores, groups=groups,which=3:4, mc.cores=2)
# plot histogram of random distances
hist(perm1$permudist, xlim=c(0,0.1),main="measured vs. random distances",
     xlab="distances")
points(perm1$dist,50,col=2,pch=19)#actual distance
text(perm1$dist,60,label=paste("actual distance\n(p=",perm1$p.value,")"))
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
