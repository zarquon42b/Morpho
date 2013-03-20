\name{rotonto}
\alias{rotonto}
\alias{rotreverse}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{rotates, translates and scales one matrix onto an other using
  Procrustes fitting
}
\description{
rotate a matrix onto an other without loosing information about the
location of the targetmatrix
and reverse this transformations using rotreverse
}
\usage{
rotonto(x, y, scale = FALSE, signref = TRUE, reflection = TRUE, weights
= NULL, centerweight = FALSE)
rotreverse(mat,rot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{k x m matrix to be rotated onto (targetmatrix)
%%     ~~Describe \code{x} here~~
}
  \item{y}{k x m matrix which will be rotated (reference matrix)
%%     ~~Describe \code{y} here~~
}
  \item{scale}{logical: scale matrix to minimize sums of squares
%%     ~~Describe \code{scale} here~~
}
  \item{signref}{logical: report if reflections were involved in the rotation
%%     ~~Describe \code{signref} here~~
} \item{mat}{matrix on which the reverse transformations have to be applied
%%     ~~Describe \code{signref} here~~
}
 \item{rot}{an object resulting from the former application of rotonto
%%     ~~Describe \code{signref} here~~
}
\item{reflection}{allow reflections.

%%     ~~Describe \code{signref} here~~
}
\item{weights}{vector of length k, containing weights for each
  landmark. 
%%     ~~Describe \code{signref} here~~
}
\item{centerweight}{logical: if weights are defined and
  centerweigths=TRUE, the matrix will be centered according to these
  weights instead of the barycenter.

%%     ~~Describe \code{signref} here~~
}

}
\value{
%%  ~Describe the value returned
  	\item{yrot }{rotated and translated matrix}
 	 \item{Y }{centred and rotated reference matrix }
	\item{X }{centred target matrix }
	\item{trans }{vector between original position of target and centered reference (during rotation process) }
	\item{transy }{vector between original position of reference and centered reference (during rotation process) }	
	\item{gamm }{rotation matrix }	
	
%% ...
}
\references{
  Lissitz, R. W., Sch6nemann, P. H., & Lingoes, J. C. (1976). A solution to the weighted Procrustes problem
in which the transformation is in agreement with the loss function. Psychometrika, 41,547-550.

%% ~put references to the literature/web site here ~
}
\author{
Stefan Schlager}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{rotmesh.onto}}
}
\examples{
library(shapes)
lims <- c(min(gorf.dat[,,1:2]),max(gorf.dat[,,1:2]))
rot <- rotonto(gorf.dat[,,1],gorf.dat[,,2]) ### rotate the second onto the first config
plot(rot$yrot,pch=19,xlim=lims,ylim=lims) ## view result
points(gorf.dat [,,2],pch=19,col=2) ## view original config
rev1 <- rotreverse(rot$yrot,rot)
points(rev1,cex=2) ### show inversion by larger circles around original configuration

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
