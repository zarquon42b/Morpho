\name{tps3d}
\Rdversion{1.1}
\alias{tps3d}

\title{
thin plate spline mapping
}
\description{
      maps a datamatrix via thin plate spline between calculated by a reference on a target configuration in 2D and 3D
}
\usage{
tps3d(M, refmat, tarmat, lambda = 0)
}

\arguments{
  \item{M}{datamatrix - e.g. the matrix information of vertices of a given surface

}
  \item{refmat}{reference matrix - e.g. landmark configuration on a surface

}
  \item{tarmat}{target matrix - e.g. landmark configuration on a target surface

  }
   \item{lambda}{integer: regularisation parameter of the TPS.

}
}

\value{returns the warped datamatrix
  

}
\references{
Bookstein FL. 1989. Principal Warps: Thin-plate splines and the decomposition of deformations. IEEE Transactions on pattern analysis and machine intelligence 11(6).
}
\author{
Stefan Schlager
}




\seealso{
\code{\link{warp_obj}}
}



