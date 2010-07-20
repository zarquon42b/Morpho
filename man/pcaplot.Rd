\name{pcaplot}
\Rdversion{1.1}
\alias{pcaplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualization of shape change
}
\description{
  visualization of the shape changes explained by Principal components
  }
\usage{
pcaplot(refshape, PCs, Scores, pcshow = c(1, 2, 3), mag = 3, color = 4, lwd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refshape}{a k x m matrix of a reference configuration - e.g. the meanshape 
%%     ~~Describe \code{refshape} here~~
}
  \item{PCs}{a matrix containing the Principal Components (e.g. the $PCs output of the procSym function).
%%     ~~Describe \code{PCs} here~~
}
  \item{Scores}{a matrix containing the PC scores (e.g. the $PCscores output of the procSym function).
%%     ~~Describe \code{Scores} here~~
}
  \item{pcshow}{a vector containing the PCscores to be visualized.
%%     ~~Describe \code{pcshow} here~~
}
  \item{mag}{a vector or an integer containing which standard deviation of which PC has to be visualized.
%%     ~~Describe \code{mag} here~~
}
  \item{color}{color of the 3d points/spheres.
%%     ~~Describe \code{color} here~~
}
  \item{lwd}{width of the lines representing the shape change.
%%     ~~Describe \code{lwd} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{procSym}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
