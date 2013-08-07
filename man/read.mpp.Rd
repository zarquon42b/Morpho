\name{read.mpp}
\alias{read.mpp}
\title{Read saved pick-points from meshlab
}
\description{
imports pick points selected with meshlab}
\usage{
read.mpp(file, info=FALSE)
}
\arguments{
  \item{file}{file to import
  }
  \item{info}{logical: if TRUE, addtional infos are returned
  }
}

\value{
  if \code{info=FALSE}:

  a matrix containing picked-points coordinates

  if \code{info=TRUE}: a list containing
  
  \item{data }{matrix containing coordinates}
  \item{info }{additional info contained in file}

}

\author{
  Stefan Schlager

}



\seealso{
\code{\link{read.pts}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
