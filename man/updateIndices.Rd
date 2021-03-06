% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updateIndices.r
\name{updateIndices}
\alias{updateIndices}
\title{update a vector of indices after removal of some referenced items}
\usage{
updateIndices(x, ignore, indexrange)
}
\arguments{
\item{x}{vector containing indices (e.g. to matrix rows)}

\item{ignore}{integer vector: remove those items from the original structure}

\item{indexrange}{maximum range of the index in the referenced item structure}
}
\description{
update a vector of indices after removal of some referenced items
}
\examples{
refItem <- matrix(1:10,5,2)
index <- c(1,3,5) # this indexes some rows of the matrix we are interested in
## now we want to ignore row 2 and 5 and want to update the index so it will still fit
indexNew <- updateIndices(index,c(2,5),indexrange=5)

## Here a more useful example:
data(boneData)
left <- c(4,6,8)
  ## determine corresponding Landmarks on the right side:
    # important: keep same order
    right <- c(3,5,7)
    pairedLM <- cbind(left,right)
## now we want to remove some landmarks and need to updated the pairedLM indices
ignore <- c(5,6)
mynewboneLM <- boneLM[-ignore,,]
pairedLMnew <- apply(pairedLM,2,updateIndices,ignore=ignore,indexrange=dim(boneLM)[1])
}
