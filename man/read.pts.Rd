\name{read.pts}
\alias{read.pts}
\title{
  reads pts files
}
\description{
  reads Landmark data exported from the software Landmark from http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\usage{
read.pts(file = "x", na = 9999)
}
\arguments{
  \item{file}{pts file 
  }
  \item{na }{specifies a value that indicates missing values
  }
}
\value{
  \item{matrix}{A Matrix containing landmark information rownames will be the names given to the landmarks in Landmark}
}


  \seealso{
    \code{\link{read.pts}}
  }
  \examples{
data(nose)
write.pts(shortnose.lm, filename="shortnose")
data <- read.pts("shortnose.pts")
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
