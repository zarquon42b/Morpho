\name{regdist.raw}
\alias{regdist.raw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{correlation between shape space and tangent space

}
\description{performs a full procrustes superimposition of landmark data and calculates the correlation between tangent and shape space.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regdist.raw(dataarray, scale = TRUE, plot = TRUE, main = "", rho =
"angle", dist.mat.out = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataarray}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
%%     ~~Describe \code{dataarray} here~~
}
\item{scale}{Logical: if TRUE, a full GPA is performed.
%%     ~~Describe \code{plot} here~~
}
  \item{plot}{Logical: whether to plot the distances between observations.
%%     ~~Describe \code{plot} here~~
}
  \item{main}{character string: Title of the plot.
%%     ~~Describe \code{main} here~~
}
  \item{rho}{chose how to calculate distances in shape space. Options: "riemdist"=Riemannian distance (function from the shapes package-takes along time to calculate), "angle"=calculates the angle between shape vectors, "sindist"=sinus of length of residual vector between shape vectors.
%%     ~~Describe \code{rho} here~~
}
  \item{dist.mat.out}{Logical: If TRUE, output will contain distance matrices.
%%     ~~Describe \code{dist.mat.out} here~~
}
%%  \item{reflect}{Logical: indicating if reflection is required for Procrustes superimposition.
%%     ~~Describe \code{reflect} here~~
%%}
}

\value{
\item{cor}{correlation coefficient between distances in shape space and tangent space}
\item{procSS }{Procrustes Sums of Squares (of full procrustes distance.}
\item{tanSS }{Tangent Sums of Squares}
\item{rhoSS }{Procrustes Sums of Squares (of angle) }
\item{euc.dist }{distance matrix of euclidean distance in Tangent space}
\item{proc.dist }{distance matrix of Procrustes distance in Shape space}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Stefan Schlager
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{regdist}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
