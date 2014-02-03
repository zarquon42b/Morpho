\name{relaxLM}
\alias{relaxLM}
\title{relax one specific 3D landmark configuration against a reference}
\usage{
relaxLM(lm, reference, SMvector, outlines = NULL, surp = NULL,
  sur.name = NULL, mesh = NULL, tol = 1e-05, deselect = FALSE,
  inc.check = TRUE, iterations = 0, fixRepro = TRUE)
}
\arguments{
  \item{lm}{k x 3 matrix containing landmark data to be
  slidden.}

  \item{reference}{k x 3 matrix containing landmark of the
  reference}

  \item{SMvector}{A vector containing the landmarks on the
  curve(s) that are allowed to slide}

  \item{outlines}{A vector (or if threre are several
  curves) a list of vectors (containing the rowindices) of
  the (Semi-)landmarks forming the curve(s) in the
  successive position on the curve - including the
  beginning and end points, that are not allowed to slide.}

  \item{surp}{A vector containing Semilandmarks positioned
  on surfaces.}

  \item{sur.name}{character: containing the filename of the
  corresponding surface.When specified, mesh has to be
  NULL.}

  \item{mesh}{triangular mesh of class "mesh3d" loaded into
  the R workspace, when specified, "sur.name" has to be
  NULL. The function \code{\link{closemeshKD}} will be used
  for reprojection onto the surface.}

  \item{tol}{numeric: Threshold for convergence in the
  sliding proces. Full Procrustes distance between actual
  result and previous iteration.}

  \item{deselect}{Logical: if TRUE, the SMvector is
  interpreted as those landmarks, that are not allowed to
  slide.}

  \item{inc.check}{Logical: if TRUE, the program stops when
  convergence criterion starts increasing and reports
  result from last iteration.}

  \item{iterations}{integer: maximum amounts the algorithm
  runs - even when 'tol' is not reached. When iterations=0,
  the algorithm runs until convergence.}

  \item{fixRepro}{logical: if \code{TRUE}, fix landmarks
  will also be projected onto the surface. If you have
  landmarks not on the surface, select
  \code{fixRepro=FALSE}}
}
\value{
returns kx3 matrix of slidden landmarks
}
\description{
relax one specific 3D landmark configuration against a
reference (e.g. a sample mean)
}
\examples{
require(rgl)
data(nose)
### relax shornose against longnose

# define fix landmarks
fix <- c(1:5,20:21)
# define surface patch by specifying row indices of matrices
# all except those defined as fix
surp <- c(1:dim(shortnose.lm)[1])[-fix]
## to reduce this example's computation time,
# we only use the right hand semi-landmarks
# (which keeps the left hand ones fix)
surp <- surp[1:316]

relax <- relaxLM(shortnose.lm[1:323, ],
         longnose.lm[1:323, ], mesh=shortnose.mesh, iterations=1,
         SMvector=fix, deselect=TRUE, surp=surp)

\dontrun{
# visualize differences red=before and green=after sliding
deformGrid3d(shortnose.lm[1:323, ], relax, ngrid=0)
# add surface
wire3d(shortnose.mesh, col="white")
}
}
\author{
Stefan Schlager
}
\references{
Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005.
Semilandmarks in Three Dimensions, in Modern Morphometrics
in Physical Anthropology. Edited by D. E. Slice, pp. 73-98.
New York: Kluwer Academic/Plenum Publishers.
}
\seealso{
\code{\link{slider3d}}
}
\keyword{~kwd1}
\keyword{~kwd2}

