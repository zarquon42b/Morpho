\name{render.meshDist}
\alias{render.meshDist}
\alias{export.meshDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  plot or save the results of meshDist 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{render}{meshDist}(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,output=FALSE,uprange=NULL,...)
\method{export}{meshDist}(x,file="default",imagedim="100x800",...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class meshDist
    %%     ~~Describe \code{x} here~~
  }
  \item{from}{numeric: minimum distance to color; default is set to 0 mm
    %%     ~~Describe \code{from} here~~
  }
  \item{to}{numeric: maximum distance to color; default is set to the maximum distance
    %%     ~~Describe \code{to} here~~
  }
  \item{steps}{integer: determines how many intermediate colors the color ramp has.
    %%     ~~Describe \code{steps} here~~
  }
  \item{ceiling}{logical: if TRUE, the next larger integer of "to" is
    used
    %%     ~~Describe \code{ceiling} here~~
  }
  \item{output}{logical: return updated colored mesh.
    %%     ~~Describe \code{ceiling} here~~
  }
  \item{uprange}{numeric between 0 and 1: restricts "to" to a quantile of
    "to", if to is NULL.
    %%     ~~Describe \code{uprange} here~~
  }
  \item{file}{character: filename for mesh and image files
    produced. E.g. "mydist" will produce the files mydist.ply and mydist.png
    %%     ~~Describe \code{file} here~~
  } 
  \item{imagedim}{character of type 100x200 where 100 determines the
    width and 200 the height of the image.
  }
  \item{\dots}{Additional parameters which will be passed to the
    methods.
  }
}
\details{
  render.meshDist renders the colored mesh and displays the color ramp.
  export.meshDist exports the colored mesh as ply file and the color
  chart as png file.
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Stefan Schlager}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{meshDist}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
