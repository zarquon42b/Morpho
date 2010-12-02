\name{patch2surf}
\alias{mc.patch2surf}
\alias{patch2surf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ maps a set of semilandmarks from a reference to the target

}
\description{maps a set of semilandmarks from a reference to the target:
1. warps the reference configuration using the "fix" landmark information in the vicinity of the target
2. find the closest point on the target mesh - either along the normal of the warped reference mesh.
or simply the closest point.


}
\usage{
patch2surf(patch, lm.ref, lm.tar, tarmesh, ray = FALSE, refmesh = NULL, tol = NULL, clean = TRUE)
mc.patch2surf(patch, lm.ref, lm.tar, tarmesh, ray = FALSE, refmesh = NULL, tol = NULL, clean = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patch}{k x 3 matrix containing the fix and the semilandmarks
%%     ~~Describe \code{patch} here~~
}
  \item{lm.ref}{matrix containing the set of fix landmark on the reference surface.
%%     ~~Describe \code{lm.ref} here~~
}
  \item{lm.tar}{matrix containing the set of fix landmark on the target surface.
%%     ~~Describe \code{lm.tar} here~~
}
  \item{tarmesh}{either a object of class mesh3d or path to a 3D mesh file (PLY,STL,OBJ) .
%%     ~~Describe \code{tarmesh} here~~
}
  \item{ray}{logical: request search for closest point on target along normals (both directions)
%%     ~~Describe \code{ray} here~~
}
  \item{refmesh}{either a object of class mesh3d or path to a 3D mesh file (PLY,STL,OBJ) . only required if ray=TRUE
%%     ~~Describe \code{refmesh} here~~
}
  \item{tol}{tolerance
%%     ~~Describe \code{tol} here~~
}
  \item{clean}{logical: clean temporary files
%%     ~~Describe \code{clean} here~~
}
}
\details{command line tools triray_project and trimesh_project from the Auxiliaries section in http://sourceforge.net/projects/morpho-rpackage/files/ installed.
}
\value{returns a k x 3 matrix of projected points.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Stefan Schlager
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
function(patch,lm.ref,lm.tar,tarmesh,ray=FALSE,refmesh=NULL,tol=NULL,clean=TRUE)
	
{ 	if (is.null(tol))
		{tol<-1e3
		}
	est0<-tps3d(patch,lm.ref,lm.tar)
	if (ray)	
		{if (is.character(refmesh))
			{refmesh_in<-file2mesh(refmesh)
			est.mesh<-mc.warp.mesh(refmesh_in,lm.ref,lm.tar,updateNormals=FALSE)
			mesh2ply(est.mesh,"est.mesh")
			patch_upd<-proj.read(est0,"est.mesh.ply",clean=FALSE)
			system(paste("triray_project out_cloud.ply",tarmesh,tol,sep=" "))
			patchpro<-ply2mesh("out_ray.ply")
			}
		
		}
	else
	{
	patchpro<-proj.read(est0,tarmesh,readnormals=F)
	}
	if (clean)
		{unlink(c("est.mesh.ply","out_cloud.ply","out_ray.ply"))
		}
	return(patchpro)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
