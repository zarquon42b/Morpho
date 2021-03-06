% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAC.r
\name{CAC}
\alias{CAC}
\title{calculate common allometric component}
\usage{
CAC(x, size, groups = NULL, log = FALSE)
}
\arguments{
\item{x}{datamatrix (e.g. with PC-scores) or 3D-array with landmark coordinates}

\item{size}{vector with Centroid sizes}

\item{groups}{grouping variable}

\item{log}{logical: use \code{log(size)}}
}
\value{
\item{CACscores}{common allometric component scores}
\item{CAC}{common allometric component}
\item{x}{(group-) centered data}
\item{sc}{CAC reprojected into original space by applying \code{CAC \%*\% x}}
\item{RSCscores}{residual shape component scores}
\item{RSC}{residual shape components}
\item{gmeans}{groupmeans}
\item{CS}{the centroid sizes (log transformed if \code{log = TRUE})}
}
\description{
calculate common allometric component
}
\examples{
data(boneData)
proc <- procSym(boneLM)
pop.sex <- name2factor(boneLM,which=3:4)
cac <- CAC(proc$rotated,proc$size,pop.sex)
plot(cac$CACscores,cac$size)#plot scores against Centroid size
cor.test(cac$CACscores,cac$size)#check for correlation
#visualize differences between large and small on the sample's consensus
\dontrun{
large <- restoreShapes(max(cac$CACscores),cac$CAC,proc$mshape)
small <- restoreShapes(min(cac$CACscores),cac$CAC,proc$mshape)
deformGrid3d(small,large,ngrid=0)
}
}
\references{
Mitteroecker P, Gunz P, Bernhard M, Schaefer K, Bookstein FL. 2004. Comparison of cranial ontogenetic trajectories among great apes and humans. Journal of Human Evolution 46(6):679-97.
}
