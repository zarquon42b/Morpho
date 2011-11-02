\name{pcaplot3d}
\Rdversion{1.1}
\alias{pcaplot3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualization of shape change
}
\description{
  visualization of the shape changes explained by Principal components
  }
\usage{
pcaplot3d.nosymproc(x, pcshow = c(1, 2, 3), mag = 3, color = 4, lwd = 1)
pcaplot3d.symproc(x, pcshow = c(1, 2, 3), mag = 3, color = 4, lwd = 1,sym=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a object derived from the function mc.procSym or procSym
    calculated on 3D coordinates.
%%     ~~Describe \code{refshape} here~~
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
