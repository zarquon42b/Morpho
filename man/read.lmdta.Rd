\name{read.lmdta}
\alias{read.lmdta}
\title{read dta files
}
\description{
reads .dta files created by the software Landmark http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\usage{
read.lmdta(file = "x")
}
\arguments{
  \item{file}{a dta file 
}
}

\value{
  \item{arr}{array containing landmarks dimnames will be Information of ID and landmark names specified in Landmark}
 \item{info }{Information extracted from the header of the dta file}
 \item{idnames }{character vector containing the names of the individuals as specified in the dta file}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
