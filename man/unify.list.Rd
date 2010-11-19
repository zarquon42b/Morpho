\name{unify.list}
\alias{mc.unify.list}
\alias{unify.list}


\title{unify meshes to gain corresponding vertices

}
\description{unifys a list of meshes using mc.unify.list
}
\usage{
	mc.unify.list(meshlist, proc, ref) ## faster mulitcore version - only for UNIX systems.
	unify.list(meshlist, proc, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meshlist}{list containing meshes to be unified, names should be identical to dimnames(proc)[[3]]
%%     ~~Describe \code{meshlist} here~~
}
  \item{proc}{array of procrustes superimposed data
%%     ~~Describe \code{proc} here~~
}
  \item{ref}{integer: specifies the position of the reference configuration in proc.

}
}
\details{
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
