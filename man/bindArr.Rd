\name{bindArr}
\alias{bindArr}
\title{concatenate multiple arrays/matrices
}
\description{
concatenate multiple 3-dimensional arrays and/or 2-dimensional matrices to one big array

}
\usage{
bindArr(..., along = 1)
}

\arguments{
  \item{\dots}{matrices and/or arrays with appropriate dimensionality to
  combine to one array.

}
  \item{along}{dimension along which to concatenate.
}
}


\value{returns array of combined matrices/arrays
}


\seealso{\code{\link{cbind}}, \code{\link{rbind}}, \code{\link{array}}
}
\examples{
A <- matrix(rnorm(18),6,3)
B <- matrix(rnorm(18),6,3)
C <- matrix(rnorm(18),6,3)

#combine to 3D-array
newArr <- bindArr(A,B,C,along=3)
#combine along first dimension
newArr2 <- bindArr(newArr,newArr,along=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
