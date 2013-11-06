\name{read.lmdta}
\alias{read.lmdta}
\title{read dta files}
\usage{
  read.lmdta(file = "x", na = 9999)
}
\arguments{
  \item{file}{a dta file}

  \item{na}{specifies a value that indicates missing
  values}
}
\value{
  \item{arr }{array containing landmarks dimnames will be
  Information of ID and landmark names specified in
  Landmark} \item{info }{Information extracted from the
  header of the dta file} \item{idnames }{character vector
  containing the names of the individuals as specified in
  the dta file}
}
\description{
  reads .dta files created by the software Landmark
  http://graphics.idav.ucdavis.edu/research/EvoMorph
}
\keyword{~kwd1}
\keyword{~kwd2}

