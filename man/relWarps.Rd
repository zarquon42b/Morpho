\name{relWarps}
\alias{relWarps}
\title{
calculate relative Warp analysis}
\description{
After Procrustes registration the data is scaled by the bending energy
 or its inverse to emphasize global/local differences when exploring a
 sample's shape. }
\usage{
relWarps(data, scale = TRUE, CSinit = TRUE, alpha = 1, tol = 1e-10, orp = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
%%     ~~Describe \code{data} here~~
}
  \item{scale}{Logical: indicating if scaling is requested
%%     ~~Describe \code{scale} here~~
}
  \item{CSinit}{Logical: if TRUE, all configurations are initially scaled to Unit Centroid Size.

}
  \item{alpha}{integer: power of the bending energy matrix. If alpha = 0 then standard Procrustes PCA is carried out. If alpha = 1 then large scale differences are emphasized, if alpha = -1 then small scale variations are emphasised.
}
  \item{tol}{tolerance for the eigenvalues of the bending energy matrix
    to be zero
}
  \item{orp}{logical: request orthogonal projection into tangent space.
%%     ~~Describe \code{orp} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{bescores }{relative warp scores}
\item{uniscores }{uniform scores}
\item{Var }{non-affine variation explained by each relative warp}
\item{mshape}{sample's conensus shape}
\item{rotated}{Procrustes superimposed data}
\item{bePCs}{Vector basis of nonaffine shape variation- relative warps}
\item{uniPCs}{Vector basis of affine shape variation - uniform component}


}
\references{Bookstein FL 1989. Principal Warps: Thin-plate splines and the decomposition of deformations. IEEE Transactions on pattern analysis and machine intelligence 11.
Bookstein FL, 1991. Morphometric tools for landmark data. Geometry and
biology. Cambridge Univ. Press, Cambridge.

Rohlf FJ,  Bookstein FL 2003. Computing the Uniform Component of Shape Variation. Systematic Biology 52:66-69.
%% ~put references to the literature/web site here ~
}
\author{Stefan Schlager
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line