\name{showPC}
\alias{showPC}
\title{
convert PCs to landmark configuration}
\description{
convert PC-scores to landmark coordinates}
\usage{
showPC(scores, PC, mshape)
}
\arguments{
  \item{scores}{vector of PC-scores
}
  \item{PC}{Principal components (eigenvectors of the covariance matrix)
    associated with 'scores'.
}
  \item{mshape}{meanshape (that was used to center the data by the PCA)
}
}
\details{Rotates and translates PC-scores derived from shape data back
  into configuration space.
}
\value{returns matrix containing landmarks
}

\author{
Stefan Schlager
}

\seealso{
\code{\link{prcomp}}, \code{\link{procSym}}
}
\examples{
library(shapes)
## generate landmarks using
##the first PC-score of the first specimen

proc <- procSym(gorf.dat)
lm <- showPC(proc$PCscores[1,1],proc$PCs[,1],proc$mshape)
plot(lm,asp=1)

##now the first 3 scores
lm2 <- showPC(proc$PCscores[1,1:3],proc$PCs[,1:3],proc$mshape)
points(lm2,col=2)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
