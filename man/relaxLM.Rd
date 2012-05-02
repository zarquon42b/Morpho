\name{relaxLM}
\alias{relaxLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
relax one specific 3D landmark configuration against a reference
}
\description{
  relax one specific 3D landmark configuration against a reference (e.g. a
sample mean)

}
\usage{
relaxLM(lm, reference, SMvector, outlines = NULL, surp = NULL, sur.name = NULL, tol = 1e-05, deselect = FALSE, inc.check = TRUE, iterations = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lm}{k x 3 matrix containing landmark data to be slidden.
%%     ~~Describe \code{lm} here~~
}
  \item{reference}{k x 3 matrix containing landmark of the reference
%%     ~~Describe \code{reference} here~~
}
  \item{SMvector}{A vector containing the landmarks on the curve(s) that are allowed to slide
%%     ~~Describe \code{SMvector} here~~
}
  \item{outlines}{A vector (or if threre are several curves) a list of vectors (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in the successive position on the curve - including the beginning and end points, that are not allowed to slide.
%%     ~~Describe \code{outlines} here~~
}
  \item{surp}{A vector containing Semilandmarks positioned on surfaces.
%%     ~~Describe \code{surp} here~~
}
  \item{sur.name}{character: containing the filename of the corresponding surface. 
%%     ~~Describe \code{sur.name} here~~
}
  \item{tol}{numeric: Threshold for convergence in the sliding
    proces. Full Procrustes distance between actual result and previous iteration.
%%     ~~Describe \code{tol} here~~
}
  \item{deselect}{Logical: if TRUE, the SMvector is interpreted as those landmarks, that are not allowed to slide.
%%     ~~Describe \code{deselect} here~~
}
  \item{inc.check}{
%%     ~~Describe \code{inc.check} here~~
}
  \item{iterations}{Logical: if TRUE, the program stops when convergence criterion starts increasing and reports result from last iteration.
%%     ~~Describe \code{iterations} here~~
}
}
\details{needs trimesh_project to be installed (http://sourceforge.net/projects/morpho-rpackage/files/)

}
\value{returns kx3 matrix of slidden landmarks
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{\code{Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. «Semilandmarks in Three Dimensions,» in Modern Morphometrics in Physical Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer Academic/Plenum Publishers.}
%% ~put references to the literature/web site here ~
}
\author{Stefan Schlager
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{slider3d}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line