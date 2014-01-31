\name{pcaplot3d}
\alias{pcaplot3d}
\alias{pcaplot3d.nosymproc}
\alias{pcaplot3d.symproc}
\title{visualization of shape variation}
\usage{
pcaplot3d(x, ...)

\method{pcaplot3d}{symproc}(x, pcshow = c(1, 2, 3), mag = 3, color = 4,
  lwd = 1, sym = TRUE, ...)

\method{pcaplot3d}{nosymproc}(x, pcshow = c(1, 2, 3), mag = 3, color = 4,
  lwd = 1, ...)
}
\arguments{
  \item{x}{a object derived from the function procSym
  calculated on 3D coordinates.}

  \item{pcshow}{a vector containing the PCscores to be
  visualized.}

  \item{mag}{a vector or an integer containing which
  standard deviation of which PC has to be visualized.}

  \item{color}{color of the 3d points/spheres.}

  \item{lwd}{width of the lines representing the shape
  change.}

  \item{sym}{logical: if TRUE the symmetric component of
  shape is displayed. Otherwise the asymmetric one.}

  \item{\dots}{Additional parameters which will be passed
  to the methods.}
}
\description{
visualization of shape change
}
\details{
visualization of the shape changes explained by Principal
components
}
\examples{
\dontrun{
data(nose)
#make a tiny sample
nosearr <- bindArr(longnose.lm, shortnose.lm, along=3)
proc <- procSym(nosearr)
pcaplot3d(proc,pcshow=1,mag=-3)#only one PC available
}
}
\seealso{
\code{\link{procSym}}
}

