\name{histGroup}
\alias{histGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot histogram for multiple groups.
}
\description{
plot a histogram for multiple groups, each group colored individually
}
\usage{
histGroup(data, groups, main = paste("Histogram of", dataname), xlab =
dataname, ylab, col = NULL, alpha = 0.5,breaks = "Sturges", legend=TRUE,legend.x=80,legend.y=80,legend.pch=15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{vector containing data.
%%     ~~Describe \code{data} here~~
}
  \item{groups}{grouping factors
%%     ~~Describe \code{groups} here~~
}
  \item{main,xlab,ylab}{these arguments to title have useful defaults here.
%%     ~~Describe \code{main} here~~
}
\item{col}{vector containing color for each group. If NULL, the function
  "rainbow" is called.
%%     ~~Describe \code{alpha} here~~
}  
 
\item{alpha}{numeric between 0 and 1. Sets the transparency of the colors
%%     ~~Describe \code{alpha} here~~
}  
 \item{breaks}{one of:
    \itemize{
      \item a vector giving the breakpoints between histogram cells,
      \item a single number giving the number of cells for the histogram,
      \item a character string naming an algorithm to compute the
      number of cells (see \sQuote{Details}),
      \item a function to compute the number of cells.
    }
    In the last three cases the number is a suggestion only.
  }
  \item{legend}{logical: if TRUE, a legend is plotted
%%     ~~Describe \code{alpha} here~~
}  
\item{legend.x}{x position of the legend from the upper left corner
%%     ~~Describe \code{alpha} here~~
}  
\item{legend.y}{y position of the legend from the upper left corners
%%     ~~Describe \code{alpha} here~~
}  
\item{legend.pch}{integer: define the symbol to visualise group colors (\code{\link{points}}
%%     ~~Describe \code{alpha} here~~
}  

}
\details{Just a wrapper for the function hist from the "graphics" package
%%  ~~ If necessary, more details than the description above ~~
}

\author{
Stefan Schlager}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{hist}}
}
\examples{
data(iris)
histGroup(iris$Petal.Length,iris$Species)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
