\name{lineplot}
\alias{lineplot}
\title{plot lines between landmarks
}
\description{
add lines connecting landmarks to visualise a sort of wireframe}
\usage{
lineplot(x, point, col = 1, lwd = 1, line_antialias = FALSE, add = TRUE)
}
\arguments{
  \item{x}{
}
  \item{point}{vector or list of vectors containing rowindices of x,
    determining which landmarks to connect.
}
  \item{col}{color of lines
}
  \item{lwd}{line width
}
  \item{line_antialias}{logical: smooth lines
}
  \item{add}{logical: add to existing plot
}
}

\author{
Stefan Schlager}
\note{
works with 2D and 3D configurations}


\seealso{
\code{\link{pcaplot3d}}
}
\examples{
library(shapes)
##2D example
plot(gorf.dat[,,1],asp=1)
lineplot(gorf.dat[,,1],point=c(1,5:2,8:6,1),col=2)

##3D example
data(nose)
points3d(shortnose.lm[1:9,])
lineplot(shortnose.lm[1:9,],point=list(c(1,3,2),c(3,4,5),c(8,6,5,7,9)),col=2)

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
