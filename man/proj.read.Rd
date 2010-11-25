\name{proj.read}
\alias{proj.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Project points onto the closest point on a mesh

}
\description{project points onto a given surface and return projected points and normals.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
proj.read(lm, mesh, readnormals = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lm}{m x 3 matrix containing 3D coordinates.
%%     ~~Describe \code{lm} here~~
}
  \item{mesh}{character: specify path to mesh file.
%%     ~~Describe \code{mesh} here~~
}
  \item{readnormals}{logical: return normals of projected points.
%%     ~~Describe \code{readnormals} here~~
}
}
\details{
file2mesh needs the command line tool trimesh_project (http://sourceforge.net/projects/morpho-rpackage/files/) installed.
}
\value{either a m x3 matrix containing projected points or a list, where
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{vb }{3 x m matrix containing projected points}
\item{normals }{3 x m matrix containing normals}
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(lm,mesh,readnormals=TRUE)
{	if (is.character(mesh))
		{proj.back(lm,mesh)
		}
	
	else 
		{mesh2ply(mesh,"dump0")
		proj.back(lm,"dump0.ply")
		unlink("dump0.ply")
		}
	
	data<-ply2mesh("out_cloud.ply",readnormals=readnormals)
	unlink("out_cloud.ply")
	return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
