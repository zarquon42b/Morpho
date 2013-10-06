\name{pls2B}
\alias{pls2B}
\title{Two-Block partial least square regression.
}
\description{
Performs a Two-Block PLS on two sets of data and assesses the
significance of each score by permutation testing}
\usage{
pls2B(x, y, tol = 1e-12, same.config=FALSE, rounds = 0, mc.cores = detectCores())
}
\arguments{
  \item{y}{array containing superimposed landmark data of the first
    block. Matrices are also allowed but the option 'same.config' will
    not work.
}
  \item{x}{array containing superimposed landmark data second block.Matrices are also allowed but the option 'same.config' will
    not work.
}
  \item{tol}{threshold for discarding singular values.
}
 \item{same.config}{logical: if \code{TRUE} each permutation includes new
   superimposition of permuted landmarks. This is necessary if both
   blocks originate from landmarks that are superimposed together.
}
  \item{rounds}{rounds of permutation testing.
}
  \item{mc.cores}{integer: determines how many cores to use for the computation. The default is autodetect. But in case, it doesn't work as expected cores can be set manually. Parallel processing is disabled on Windows due to occasional errors. 

}
}
\details{The Two-Block PLS tries to find those linear combinations in
  each block maximising the covariance between blocks. The significance
  of each linear combination is assessed by comparing the singular value
  to those obtained from permuted blocks. If both blocks contain
  landmarks superimposed TOGETHER, the option \code{same.config=TRUE}
  requests superimposition of the permuted configurations (i.e. where
  the the landmarks of block \code{x} are replaced by corresponding landmarks
  of other specimen.
}
\value{

  \item{svd }{singular value decomposition (see \code{\link{svd}}) of the 'common' covariance block.}
  \item{Xscores }{PLS-scores of x}
  \item{Yscores }{PLS-scores of y}
   \item{CoVar }{Dataframe containing singular values, explained
  covariation, correlation coeffictient between PLS-scores and p-values.}
}
\references{Rohlf FJ, Corti M. 2000. Use of two-block partial least-squares to study covariation in shape. Systematic Biology 49:740-753.


}
\author{
Stefan Schlager}



\seealso{
 \code{\link{svd}}
}
\examples{
library(shapes)
### very arbitrary test:
### check if first 4 landmarks covaries with the second 4
proc <- procSym(gorf.dat)
## we do only 50 rounds to minimize computation time
pls1 <- pls2B(proc$rotated[1:4,,],proc$rotated[5:8,,],
              same.config=TRUE,rounds=50,mc.cores=2)
pls1$CoVar
layout(matrix(1:4,2,2,byrow=TRUE))
for(i in 1:4)
 plot(pls1$Xscores[,i]~pls1$Yscores[,i])



}
\keyword{dynamic}
