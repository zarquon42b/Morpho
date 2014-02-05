\name{covDist}
\alias{covDist}
\alias{covPCA}
\title{calculate distances and PC-coordinates of covariance matrices
\code{covDist} calculates the Distance between covariance matrices while \code{covPCA} uses a MDS (multidimensional scaling) approach to obtain PC-coordinates
from a  distance matrix derived from multiple groups. P-values for pairwise
distances can be computed by permuting group membership and comparing actual
distances to those obtained from random resampling.}
\usage{
covDist(s1, s2)

covPCA(data, groups, rounds = 1000)
}
\arguments{
  \item{s1}{m x m covariance matrix}

  \item{s2}{m x m covariance matrix}

  \item{data}{matrix containing data with one row per
  observation}

  \item{groups}{factor: group assignment for each specimen}

  \item{rounds}{integer: rounds to run permutation of
  distances by randomly assigning group membership}
}
\value{
\code{covDist} returns the distance between s1 and s2

\code{covPCA} returns a list containing:

if \code{scores = TRUE} \item{PCscores }{PCscores}
\item{eigen}{eigen decomposition of the centered inner
product} if \code{rounds > 0} \item{dist }{distance matrix}
\item{p.matrix }{p-values for pairwise distances from
permutation testing}
}
\description{
calculate distances and PC-coordinates of covariance
matrices
}
\details{
calculate PC-coordinates of covariance matrices by using
the Riemannian metric in their respective space.
}
\examples{
cpca <- covPCA(iris[,1:4],iris[,5])
cpca$p.matrix #show pairwise p-values for equal covariance matrices
\dontrun{
require(car)
sp(cpca$PCscores[,1],cpca$PCscores[,2],groups=levels(iris[,5]),
   smooth=FALSE,xlim=range(cpca$PCscores),ylim=range(cpca$PCscores))

data(boneData)
proc <- procSym(boneLM)
pop <- name2factor(boneLM, which=3)
## compare covariance matrices for PCscores of Procrustes fitted data
cpca1 <- covPCA(proc$PCscores, groups=pop, rounds = 1000)
## view p-values:
cpca1$p.matrix # differences between covariance matrices
# are significant
## visualize covariance ellipses of first 5 PCs of shape
spm(proc$PCscores[,1:5], groups=pop, smooth=FALSE,ellipse=TRUE, by.groups=TRUE)
## covariance seems to differ between 1st and 5th PC
## for demonstration purposes, try only first 4 PCs
cpca2 <- covPCA(proc$PCscores[,1:4], groups=pop, rounds = 1000)
## view p-values:
cpca2$p.matrix # significance is gone
}
}
\author{
Stefan Schlager
}
\references{
Mitteroecker P, Bookstein F. 2009. The ontogenetic
trajectory of the phenotypic covariance matrix, with
examples from craniofacial shape in rats and humans.
Evolution 63:727-737.

Hastie T, Tibshirani R, Friedman JJH.  2013. The elements
of statistical learning. Springer New York.
}
\seealso{
\code{\link{prcomp}}
}
\keyword{~kwd1}
\keyword{~kwd2}

