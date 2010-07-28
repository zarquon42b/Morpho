\name{CVA}
\Rdversion{1.1}
\alias{CVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Canonical Variate Analysis
}
\description{
performs a Canonical Variate Analysis.
}
\usage{
CVA(dataarray, groups, weighting = TRUE, tolinv = 1e-10, plot = TRUE,rounds=10000,cv=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataarray}{Either a k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size. Or alternatively a n x m Matrix where n is the numeber of observations and m the number of variables (this can be PC scores for example)
}
  \item{groups}{a chracter/factor vector or a list of numeric vectors containing either the positions in an array or the rows of the matrix of the observations belonging to the same group)
%%     ~~Describe \code{groups} here~~
}
  \item{weighting}{Logical: Determines whether the between group covariance matrix is to be weighted according to group size.
%%     ~~Describe \code{weighting} here~~
}
  \item{tolinv}{Threshold for the eigenvalues of the pooled within-group-covariance matrix to be taken as zero - for calculating the general inverse of the pooled withing groups covariance matrix.
%%     ~~Describe \code{tolinv} here~~
}
  \item{plot}{Logical: determins whether in the two-sample case a histogramm ist to be plotted.
%%     ~~Describe \code{plot} here~~
}
 \item{rounds}{integer: number of permutations if a permutation test of the mahalanobis distance between group means is requested.If rounds = 0, no test is performed.
%%     ~~Describe \code{plot} here~~
}
 \item{cv}{logical: requests a Jackknife Crossvalidation.
%%     ~~Describe \code{plot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use      CV = CV, CVscores = CVscores, Grandm = Grandm,    groupmeans = groupmeans, Var = Var, CVvis = CVvis, B = B
        
\item{CV }{A matrix containing the Canonical Variates}
%%  \item{comp2 }{Description of 'comp2'}
\item{CVscores }{A matrix containing the individual Canonical Variate scores}
\item{Grandm}{a vector or a matrix containing the Grand Mean (depending if the input is an array or a matrix)}
\item{groupmeans }{a matrix or an array containing the group means (depending if the input is an array or a matrix)}
\item{Var }{Variance explained by the Canonical Variates}
\item{CVvis }{Canonical Variates projected back into the original space - to be used for visualization purposes, for details see example below.}
\item{Dist }{Mahalanobis Distances between group means - if requested tested by permutation test if the input is an array it is assumed to be superimposed Landmark Data and Procrustes Distance will be calculated.}
\item{CVcv }{A matrix containing crossvalidated CV scores}
}
\references{Cambell, N. A. & Atchley, W. R.. 1981 The Geometry of Canonical Variate Analysis: Syst. Zool., 30(3), 268-280.
Klingenberg, C. P. & Monteiro, L. R. 2005 Distances and directions in multidimensional shape spaces: implications for morphometric applications. Systematic Biology 54, 678-688.
%% ~put references to the literature/web site here ~
}
\author{
Stefan Schlager
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
alldat<-procSym(abind(gorf.dat,gorm.dat))   # perform procrustes fit on raw data
groups<-list(c(1:30),c(31:59))              # create list with group affinities
cvall<-CVA(alldat$orpdata,groups,rounds=1000)     # perform CVA and test Mahalanobis distance between groups with permutation test by 1000 rounds)

### visualize the effects of the variates:
dims<-dim(alldat$mshape)## get dims of mean
plot(alldat$mshape+4*matrix(cvall$CVvis[,1],dims[1],dims [2]),col=2,asp=1) ## plot shape changes associated with Canonical Variate 1 (as there is only one in this example) and a score of 4, color = red
points(alldat$mshape+-5*matrix(cvall$CVvis[,1],dims[1],dims [2]),col=5) ## plot shape changes associated with Canonical Variate[,1] (as there is only one in this example) and a score of -5, color = blue

points(alldat$mshape,col=3) ## add meanshape to plot, color = green



 }
