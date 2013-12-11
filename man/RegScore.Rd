\name{RegScore}
\alias{RegScore}
\title{calulate regression scores for linear model}
\usage{
  RegScore(model, x = NULL)
}
\arguments{
  \item{model}{linear model}

  \item{x}{optional: matrix containing new data to be
  projected onto the regression lines.}
}
\value{
  returns a n x m matrix containing the regression scores
  for each specimen.
}
\description{
  calulate regression scores for linear model as specified
  in Drake & Klingenberg(2008)
}
\details{
  the data is orthogonally projected onto the regression
  lines associated with each factor.
}
\examples{
model <- lm(as.matrix(iris[,1:2]) ~ iris[,3])
rs <- RegScore(model)
plot(rs,iris[,4])
}
\references{
  Drake, AG. & Klingenberg, CP. The pace of morphological
  change: historical transformation of skull shape in St
  Bernard dogs. Proceedings of the Royal Society B:
  Biological Sciences, The Royal Society, 2008, 275, 71-76.
}

