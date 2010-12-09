\name{mc.unify.list}
\alias{mc.unify.list}



\title{unify meshes to gain corresponding vertices

}
\description{unifys a list of meshes using mc.unify.list
}
\usage{
	mc.unify.list(lm.data,ref,surpath="sur",ext=".ply",files=NULL,ray=FALSE,tol=NULL,rawout=NULL,uniout=NULL,keep=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lm.data}{array containing landmark data
%%     ~~Describe \code{meshlist} here~~
}
  \item{ref}{integer: specifies the position of the reference configuration in lm.data.
%%     ~~Describe \code{proc} here~~
}
  \item{surpath}{character: path to folder, that contains mesh files.

}
  \item{ext}{character: extension of files (only necessary if "meshes" is not specified).

}
 \item{files}{character vector: vector of filenames in the same order as corresponding landmark configuration in lm.data.

}
\item{ray}{logical: requests search of next vertex along normals
%%     ~~Describe \code{mesh2} here~~
}
\item{tol}{numeric: threshold for maximal distance along normals if ray=TRUE (in both directions of the normal), if nothing is hit within this threshold, the closest point on the target will be used.
%%     ~~Describe \code{mesh2} here~~
}
\item{rawout}{character: if specified, the unified mesh in the original coordinate system will be saved in the specified folder using the original name but with the file extension ".ply".
%%     ~~Describe \code{mesh2} here~~
}
\item{rawout}{character: if specified, the unified mesh in the procrustes coordinate system will be saved in the specified folder using the original name but with the file extension ".ply".
%%     ~~Describe \code{mesh2} here~~
}
\item{keep}{logical: if FALSE, the meshes will only be written to disk, and output will be an empty list.
%%     ~~Describe \code{mesh2} here~~
}
}
\details{if "files" is not specified, it is necessary that  dimnames(lm.data[[3]]) are identical to the name of the surface files (without the file extension).
E.g.: if observation nr. 3 is called obs3, the corresponding mesh has to be named "obs3.ply".
%%  ~~ If necessary, more details than the description above ~~
}
\value{a list containing unified meshes (of class "mesh3d")
%%  ~Describe the value returned
%%  If it is a LIST, use
%%\item{mesh }{Description of 'comp1'}
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(meshlist,proc,ref)
{	
	n<-dim(proc)[[3]]
	uni.list<-list()
	nam<-names(meshlist)
	namref<-dimnames(proc)[[3]][ref]
	jref<-grep(namref,nam)
	njref<-names(meshlist)[[jref]]
	
	
	for (i in 1:n)
	{	nm0<-dimnames(proc)[[3]][i]
		j<-grep(nm0,nam)
		nj<-names(meshlist)[[j]]
		if (i == ref)
			{uni.list[[i]]<-meshlist[[j]]
			
			}
		else
			{uni.list[[i]]<-mc.unify.mesh(meshlist[[jref]],meshlist[[j]],proc[,,ref],proc[,,i])$unimesh
			}
		
		
		names(uni.list)[i]<-nj
		#print(j)
	}
		
		return(uni.list)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
