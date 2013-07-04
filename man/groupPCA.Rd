\name{groupPCA}
\alias{groupPCA}
\title{
Perform PCA based of the group means' covariance matrix}
\description{
Calculate covariance matrix of the groupmeans and project all
observations into the eigenspace of this covariance matrix. This
displays a low dimensional between group structure of a high dimensional problem.
}
\usage{
groupPCA(dataarray, groups, rounds = 10000, tol = 1e-10, cv=TRUE,
mc.cores = detectCores(),weighting=TRUE)
}
\arguments{
  \item{dataarray}{Either a k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size. Or alternatively a n x m Matrix where n is the numeber of observations and m the number of variables (this can be PC scores for example)
}
  \item{groups}{a chracter/factor vector or a list of numeric vectors containing either the positions in an array or the rows of the matrix of the observations belonging to the same group)
}
  \item{rounds}{integer: number of permutations if a permutation test of the euclidean distance between group means is requested.If rounds = 0, no test is performed.
}
  \item{tol}{threshold to ignore eigenvalues of the covariance matrix.}
    
    \item{cv}{logical: requests leaving-one-out crossvalidation
    
}

  \item{mc.cores}{integer: how many cores of the Computer are allowed to be
    used. Default is use autodetection by using detectCores()
    from the parallel package. Parallel processing is disabled on
    Windows due to occasional errors.
}
 \item{weighting}{logical:weight between groups covariance matrix
    according to group sizes.
}
}

\value{
  \item{eigenvalues }{Non-zero eigenvalues of the groupmean covariance matrix.}
  \item{groupPCs }{PC-axes - i.e. eigenvectors of the groupmean covariance matrix. }
  \item{Variance}{table displaying the variance explained by eache eigenvalue. }
  \item{Scores}{Scores of all observation in the PC-space}
  \item{probs}{p-values of pairwise groupdifferences - based on
    permuation testing }
  \item{groupdists}{Euclidean distances between groups' averages}
  \item{groupmeans}{Groupmeans}
  \item{Grandmean}{Grand mean}
  \item{CV}{Cross-validated scores}



}
\references{Mitteroecker P,  Bookstein F 2011. Linear Discrimination,
  Ordination, and the Visualization of Selection Gradients in Modern
  Morphometrics. Evolutionary Biology 38:100-114.

  
Boulesteix, A. L. 2005: A note on between-group PCA, International Journal of Pure and Applied Mathematics 19, 359-366.

  
}
\author{
Stefan Schlager}



\seealso{
 \code{\link{CVA}}, ~~~
}
\examples{
require(car)
data(iris)
vari=iris[,1:4]
facto=iris[,5]
pca.1=groupPCA(vari,groups=facto,rounds=1000,mc.cores=1)

### plot scores
scatterplotMatrix(~pca.1$Scores | facto, ellipse=TRUE,
        by.groups=TRUE,var.labels=c("PC1","PC2","PC3"))

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
