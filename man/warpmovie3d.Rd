\name{warpmovie3d}
\alias{warpmovie3d.mesh3d}
\alias{warpmovie3d.matrix}
\alias{warpmovie3d}
\alias{warpmovie2d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a sequence of images showing predefined steps of warping two 
meshes or landmark configurations (2D and 3D) into each other}
\description{given two landmark configurations or two meshes with the same amount of vertices and faces
  (e.g a mesh and its warped counterpart), the starting configuration/mesh will be
  subsequently transformed into the final configuration/mesh by splitting the
  differences into a predefined set of steps.
}
\usage{
\method{warpmovie3d}{mesh3d}(x, y, n, col = "green", palindrome = FALSE,
folder = NULL, movie = "warpmovie", add = FALSE, close = TRUE,countbegin=0,ask=TRUE,...)
\method{warpmovie3d}{matrix}(x, y, n, col = "green", palindrome = FALSE,
folder = NULL, movie = "warpmovie", add = FALSE, close=TRUE,
countbegin=0, ask=TRUE, radius=NULL, links=NULL, lwd=1,...)

warpmovie2d(x, y, n, col = "green", palindrome = FALSE, folder =
NULL,movie="warpmovie", links=NULL, lwd=1, imagedim = "800x800",
par=list(xaxt="n", yaxt="n", bty="n"),...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{mesh to start with (object of class mesh3d)
%%     ~~Describe \code{x} here~~
}
  \item{y}{resulting mesh (object of class mesh3d), having the same
    amount of vertices and faces than the starting mesh
%%     ~~Describe \code{y} here~~
}
  \item{n}{integer: amount of intermediate steps.
%%     ~~Describe \code{n} here~~
}
  \item{col}{color of the mesh
%%     ~~Describe \code{col} here~~
}
  \item{palindrome}{logical: if TRUE, the procedure will go forth and back.
%%     ~~Describe \code{palindrome} here~~
}
  \item{folder}{character: output folder for created images (optional)
%%     ~~Describe \code{folder} here~~
}
  \item{movie}{character: name of the output files
%%     ~~Describe \code{movie} here~~
}
\item{add}{logical: if TRUE, the movie will be added to the focussed rgl-windows.
%%     ~~Describe \code{add} here~~
} 
  \item{close}{logical: if TRUE, the rgl window will be closed when finished.
    width and 200 the height of the image.
  }
  \item{countbegin}{integer: number to start image sequence. 
    
  }
  \item{ask}{logical: if TRUE, the viewpoint can be selected manually.
    
  }
\item{radius}{numeric: define size of spheres (overides atuomatic size estimation).
%%     ~~Describe \code{add} here~~
} 
\item{links}{vector or list of vectors containing wireframe information
  to connect landmarks (optional).
%%     ~~Describe \code{add} here~~
}
\item{lwd}{numeric: controls width of lines defined by "links".
  %%     ~~Describe \code{add} here~~
  
}
 \item{imagedim}{character of pattern "100x200" where 100 determines the
    width and 200 the height of the image.
  }
 
   \item{par}{list of graphial parameters: details can be found here: \code{\link{par}}.
  }
\item{\dots}{additional arguments passed to \code{\link{shade3d}} (3D)
 or \code{\link{points}} (2D).
%%     ~~Describe \code{add} here~~
}
}
\details{A series of png files will be saved to disk. These can be
  joined to animated gifs by external programs such as imagemagick or
  used to create animations in PDFs in a latex environment (e.g. latex package: aninmate).
}

\author{
Stefan Schlager}

\seealso{
    \code{\link{ply2mesh},\link{file2mesh},\link{mesh2ply},\link{warp.mesh}}
  }

  \examples{

###3D example
data(nose)##load data
##warp a mesh onto another landmark configuration:
warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
\dontrun{
warpmovie3d(shortnose.mesh,warpnose.long,n=15)## create 15 images.
}
### restrict to landmarks
\dontrun{
warpmovie3d(shortnose.lm,longnose.lm,n=15,movie="matrixmovie")## create 15 images.
}
### the images are now stored in your current working directory and can
### be concatenated to a gif using an external program such as
### imagemagick.

### 2D example
library(shapes)
bb <- procSym(gorf.dat)
### morph superimposed first specimen onto sample mean
warpmovie2d(bb$rotated[,,1],bb$mshape,n=20,links=c(1,5,4:2,8:6,1),imagedim="600x400")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
