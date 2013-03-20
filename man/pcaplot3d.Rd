\name{pcaplot3d}
\Rdversion{1.1}
\alias{pcaplot3d}
\alias{pcaplot3d.nosymproc}
\alias{pcaplot3d.symproc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visualization of shape change
}
\description{
  visualization of the shape changes explained by Principal components
  }
  \usage{
pcaplot3d(x,...)
\method{pcaplot3d}{symproc}(x, pcshow = c(1, 2, 3), mag = 3, color = 4, lwd = 1,sym=TRUE,...)
\method{pcaplot3d}{nosymproc}(x, pcshow = c(1, 2, 3), mag = 3, color = 4, lwd = 1,...)

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
 \item{sym}{logical: if TRUE the symmetric component of shape is
  displayed. Otherwise the asymmetric one.
%%     ~~Describe \code{lwd} here~~
}
\item{\dots}{Additional parameters which will be passed to the
  methods.
}
}


\seealso{\code{\link{procSym}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
