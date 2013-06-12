\name{covDist}
\alias{covPCA}
\alias{covDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  calculate distances and PC-coordinates of covariance matrices
}
\description{
  calculate PC-coordinates of covariance matrices by using the Riemannian
  metric in their respective space.
}
\usage{
covDist(s1, s2)
covPCA(data,groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s1}{m x m covariance matrix 
    %%     ~~Describe \code{s1} here~~
  }
  \item{s2}{m x m covariance matrix 
    %%     ~~Describe \code{s2} here~~
  }
  \item{data}{matrix containing data with one row per observation.
    %%     ~~Describe \code{s2} here~~
  }
  \item{groups}{factor: group assignment for each specimen.
    %%     ~~Describe \code{s2} here~~
  }

}
\details{
  covDist calculates the Distance between covariance matrices. And covPCA
  uses a MDS (multidimensional scaling) approach to obtain PC-coordinates
  from a distance matrix derived from multiple groups.
}
\value{
  covDist returns the distance between s1 and s2, while covPCA returns a
  list containing:
  
  \item{PCscores }{PCscores
  }
  \item{dist }{distance matrix
  }
  \item{eigen }{eigen decomposition of the centered inner product
  }

}
\references{Mitteroecker P, Bookstein F. 2009. The ontogenetic
  trajectory of the phenotypic covariance matrix, with examples from
  craniofacial shape in rats and humans. Evolution 63:727-737.

  Hastie T, Tibshirani R, Friedman JJH.  2013. The elements of statistical learning. Springer New York, .
}
\author{Stefan Schlager
 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{prcomp}}
}
\examples{
require(car)
cpca <- covPCA(iris[,1:4],iris[,5])
sp(cpca$PCscores[,1],cpca$PCscores[,2],groups=levels(iris[,5]),smooth=FALSE,xlim=range(cpca$PCscores),ylim=range(cpca$PCscores))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
