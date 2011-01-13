\name{mc.fixLMtps}
\alias{mc.fixLMtps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{estimate missing landmarks
}
\description{if landmarks are missing they will be replaced by a warping of a population mean or a weighted estimate onto the deficient configuration using a Thin-Plate-Spline based on the existing landmarks
}
\usage{
mc.fixLMtps(data,comp=3,weight=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{array containing landmark data
%%     ~~Describe \code{data} here~~
}
  \item{comp}{integer: select how many of the closest observations are to be taken to calculate an initial estimate.
%%     ~~Describe \code{data} here~~
}
  \item{weight}{logical: requests the calculation of an estimate based on the procrustes distance.
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{out }{array containing all data, including fixed configurations  - same order as input}
\item{mshape }{meanshape  - calculated from complete datasets.}
\item{checklist }{list containing information about missing landmarks}
\item{check }{vector containing position of observations in data where at least one missing coordinate was found }
%% ...
}
\references{
%
}
\author{Stefan Schlager
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
function(data)
{		n<-dim(data)[3]
		k<-dim(data)[1]
		m<-dim(data)[2]
		checklist<-list()
		checkvec<-rep(0,n)
		out<-data
	### check for missing landmarks ###
		for (i in 1:n)
			{count<-0
			found<-FALSE
			checklist[[i]]<-NA
			for (j in 1:k)
				{
				if (NA\%in\%data[j,,i])
					{count<-count+1
					checklist[[i]][count]<-j
					checkvec[i]<-1
					}
				}
			}
		#### calc mean of complete configs ###
		check<-which(checkvec==1)
		data.c<-data[,,-check]
		#print(data.c)
		mean0<-mc.procGPA(data.c)$mshape
		
		for (i in 1:length(check))
			{miss<-checklist[[check[i]]]
			out[,,check[i]]<-mc.tps3d(mean0,mean0[-miss,],data[-miss,,check[i]])
			}
		return(list(out=out,mshape=mean0,checklist=checklist,check=check))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
