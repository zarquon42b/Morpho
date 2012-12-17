\name{render.meshDist}
\alias{render.meshDist}
\alias{render.matrixDist}
\alias{export.meshDist}
\alias{render}
\alias{export}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  plot or save the results of meshDist 
%%  ~~function to do ... ~~
}
\description{  Visualise or save the results of meshDist to disk. 
}
\usage{
\method{render}{meshDist}(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,displace = FALSE, shade = TRUE,sign=NULL,...)

\method{render}{matrixDist}(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,type=c("s","p"),radius=NULL,displace = FALSE, sign =
                 NULL,...)

\method{export}{meshDist}(x,file="default",imagedim="100x800",...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class meshDist
  }
  \item{from}{numeric: minimum distance to color; default is set to 0 mm
  }
  \item{to}{numeric: maximum distance to color; default is set to the maximum distance
  }
  \item{steps}{integer: determines how many intermediate colors the color ramp has.
  }
  \item{ceiling}{logical: if TRUE, the next larger integer of "to" is
    used
  }
  \item{uprange}{numeric between 0 and 1: restricts "to" to a quantile of
    "to", if to is NULL.
  }
  \item{tol}{numeric: threshold to color distances within this threshold
    green.
  }
  \item{displace}{logical: if TRUE, displacement vectors between original
    and closest points are drawn colored according to the distance.
  }
  \item{shade}{logical: if FALSE, the rendering of the colored surface
    will be supressed.
  }
  \item{sign}{logical: request signed distances to be visualised.
  }   
  \item{file}{character: filename for mesh and image files
    produced. E.g. "mydist" will produce the files mydist.ply and mydist.png
  } 
  \item{imagedim}{character of pattern "100x200" where 100 determines the
    width and 200 the height of the image.
  }
  \item{type}{character: "s" shows coordinates as spheres, while "p" shows
    3D dots.
  }
  \item{radius}{determines size of spheres; if not specified, optimal
    radius size will be estimated by centroid size of the configuration.
  }

  \item{\dots}{for render.meshDist: additional arguments passed to  \code{\link{shade3d}}. See  \code{\link{rgl.material}} for details.
  }
}
\details{
  render.meshDist renders the colored mesh and displays the color ramp
  and returns an object of class "meshDist".
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
 \code{\link{meshDist}}, \code{\link{shade3d}}

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
