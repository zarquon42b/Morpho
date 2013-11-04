\name{covDist}
\alias{covPCA}
\alias{covDist}
\title{
  calculate distances and PC-coordinates of covariance matrices
}
\description{
  calculate PC-coordinates of covariance matrices by using the Riemannian
  metric in their respective space.
}
\usage{
covDist(s1, s2)
covPCA(data, groups, scores=TRUE, rounds=0, mc.cores=detectCores())
}
\arguments{
  \item{s1}{m x m covariance matrix 
    
  }
  \item{s2}{m x m covariance matrix 
    
  }
  \item{data}{matrix containing data with one row per observation.
    
  }
  \item{groups}{factor: group assignment for each specimen.
    
  }
  \item{scores}{logical: PCscores are calculated from distance matrix.
    
  }
  \item{rounds}{integer: rounds to run permutation of distances by
  randomly assigning group membership.
    
  }
  \item{mc.cores}{integer: how many CPU-cores shall be used in
  permutation testing (not available on Windows)
    
  }

}
\details{
  covDist calculates the Distance between covariance matrices. And covPCA
  uses a MDS (multidimensional scaling) approach to obtain PC-coordinates
  from a distance matrix derived from multiple groups. P-values for
  pairwise distances can be computed by permuting group membership and
  comparing actual distances to those obtained from random resampling.
}
\value{
  \code{covDist} returns the distance between s1 and s2, while \code{covPCA} returns a
  list containing:
   \item{dist }{distance matrix
   }

   if scores = TRUE 
  \item{PCscores }{PCscores
  }
 
  \item{eigen }{eigen decomposition of the centered inner product
  }

  if rounds > 0 
  \item{p.matrix }{p-values for pairwise distances from permutation testing
  }
}
\references{Mitteroecker P, Bookstein F. 2009. The ontogenetic
  trajectory of the phenotypic covariance matrix, with examples from
  craniofacial shape in rats and humans. Evolution 63:727-737.

  Hastie T, Tibshirani R, Friedman JJH.  2013. The elements of statistical learning. Springer New York.
}
\author{Stefan Schlager
 
}



\seealso{
  \code{\link{prcomp}}
}
\examples{
require(car)
cpca <- covPCA(iris[,1:4],iris[,5])
\dontrun{
sp(cpca$PCscores[,1],cpca$PCscores[,2],groups=levels(iris[,5]),
   smooth=FALSE,xlim=range(cpca$PCscores),ylim=range(cpca$PCscores))

data(boneData)
proc <- procSym(boneLM)
pop <- name2factor(boneLM, which=3)
## compare covariance matrices for PCscores of Procrustes fitted data
cpca1 <- covPCA(proc$PCscores, groups=pop, rounds = 100, mc.cores=2)
## view p-values:
cpca1$p.matrix # differences between covariance matrices
# are significant
## visualize covariance ellipses of first 5 PCs of shape
spm(proc$PCscores[,1:5], groups=pop, smooth=FALSE,ellipse=TRUE, by.groups=TRUE)
## covariance seems to differ between 1st and 5th PC
## for demonstration purposes, try only first 4 PCs
cpca2 <- covPCA(proc$PCscores[,1:4], groups=pop, rounds = 100,
mc.cores=2)
## view p-values:
cpca2$p.matrix # significance is gone
}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
