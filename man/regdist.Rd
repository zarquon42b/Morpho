\name{regdist}
\alias{regdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{correlation between shape space and tangent space

}
\description{calculates and visualizes correlation between shape space and tangent space to see if shape variation is small enough that shape can be analyzed in Tangent space.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
regdist(proc, plot = TRUE, main = "", rho = "angle", dist.mat.out = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{proc}{the outcome of the functions procSym or procGPA (of the shapes package).
%%     ~~Describe \code{proc} here~~
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
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use     (list(cor = correlation, procSS = procdis, tanSS = eudis,   rhoSS = procdis2, euc.dist = em, proc.dist = procvec))
            
\item{cor}{correlation coefficient between distances in shape space and tangent space}
\item{procSS }{Procrustes Sums of Squares (of full procrustes distance.}
\item{tanSS }{Tangent Sums of Squares}
\item{rhoSS }{Procrustes Sums of Squares (of angle) }
\item{euc.dist }{distance matrix of euclidean distance in Tangent space}
\item{proc.dist }{distance matrix of Procrustes distance in Shape space}
%% ...

}

\author{
Stefan Schlager
}


\seealso{
\code{\link{regdist.raw}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
