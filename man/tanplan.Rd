\name{tanplan}
\alias{tanplan}
\alias{crossp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate the orthogonal complement of a 3D-vector}
\description{
calculate the orthogonal complement of a 3D-vector or the
3D-crossproduct, finding an orthogonal vector to a plane in 3D.}

\usage{
tanplan(x)
crossp(x,y)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of length 3.
  }
  \item{y}{vector of length 3.
  }
}

\value{
  tanplan:
  
  \item{y }{vector orthogonal to x}
  \item{z }{vector orthogonal to x and y}

  crossp: returns a vector of length 3.
}

\author{
Stefan Schlager}
