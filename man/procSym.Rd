\name{mc.procSym}
\alias{mc.procSym}
\alias{procSym}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Procrustes registration
  
}
\description{
    \code{procSym} performs Procrustes superimposition including sliding of semilandmarks in 3d
}
\usage{
mc.procSym(dataarray,pairedLM=NULL,SMvector=NULL,outlines=NULL,orp=TRUE,tol=1e-05,CSinit=TRUE,deselect=FALSE,recursive=TRUE,iterations=0,scale=TRUE,sizeshape=FALSE,initproc=FALSE,use.lm=NULL,center.part=TRUE)
procSym(dataarray,pairedLM=NULL,SMvector=NULL,outlines=NULL,orp=TRUE,tol=1e-05,CSinit=TRUE,deselect=FALSE,recursive=TRUE,iterations=0,scale=TRUE,sizeshape=FALSE,initproc=FALSE,use.lm=NULL,center.part=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataarray}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size.
%%     ~~Describe \code{dataarray} here~~
}
  \item{pairedLM}{A X x 2 matrix containing the paired LM. E.g. the left column contains the lefthand landmarks, 
                  while the right side contains the corresponding right hand landmarks.
%%     ~~Describe \code{pairedLM} here~~
}
  \item{SMvector}{A vector containing the landmarks on the curve(s) that are allowed to slide
%%     ~~Describe \code{SMvector} here~~
}
  \item{outlines}{A vector (or if threre are several curves) a list of vectors (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in the successive position on the curve - including the beginning and end points, that are not allowed to slide.
%%     ~~Describe \code{outlines} here~~
}
  \item{orp}{Logical: if TRUE, an orthogonal projection into shape space is performed.
%%     ~~Describe \code{orp} here~~
}
  \item{tol}{numeric: Threshold for convergence in the sliding process
%%     ~~Describe \code{tol} here~~
}
  \item{CSinit}{Logical: if TRUE, all configurations are initially scaled to Unit Centroid Size.
%%     ~~Describe \code{CSinit} here~~
}
  \item{deselect}{Logical: if TRUE, the SMvector is interpreted as those landmarks, that are not allowed to slide.
%%     ~~Describe \code{deselect} here~~
}
  \item{recursive}{Logical: if TRUE, during the iterations of the sliding process, the outcome of the previous iteration will be used. 
                    Otherwise the original configuration will be used in all iterations. 
%%     ~~Describe \code{recursive} here~~
}
  \item{iterations}{integer: select manually how many iterations will be performed during the sliding process 
                    (usefull, when there is very slow convergence). 
                              0 means iteration until convergence.
%%     ~~Describe \code{iterations} here~~
}
  \item{scale}{Logical: indicating if scaling is requested
%%     ~~Describe \code{scale} here~~
}
  \item{sizeshape}{Logical: if TRUE, a log transformed variable of Centroid Size will be added to the shapedata as first variable before performing the PCA. 
%%     ~~Describe \code{reflect} here~~
}

  \item{initproc}{Logical: indicating if the first Relaxation step is performed against the mean of an initial Procrustes superimposition. 
                  Symmetrical configurations will be relaxed against a perfectly symmetrical mean.
%%     ~~Describe \code{initproc} here~~
}
 \item{use.lm}{vector of integers to define a subset of landmarks to be used for Procrustes registration-

 }
  \item{center.part}{Logical: if TRUE, the data superimposed by the
  subset defined by use.lm will be centered according to the centroid of
  the complete configuration. Otherwise orp will be set to FALSE to
  avoid erroneous projection into tangent space.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{size }{a vector containing the Centroid Size of the configurations}
  
  \item{rotated }{k x m x n array of the rotated configurations}
  
  \item{Sym }{k x m x n array of the Symmetrical component - only available for the "Symmetry"-Option}
  
  \item{Asym }{k x m x n array of the Asymmetrical component - only available for the "Symmetry"-Option}
  \item{asymmean }{k x m matrix of mean asymmetric deviation from symmetric mean}
  \item{mshape }{sample meanshape}
  \item{symmean }{meanshape of symmetrized configurations}
  \item{tan }{if orp=TRUE: Residuals in tangentspace else, Procrustes residuals - only available without the "Symmetrie"-Option}
  \item{PCs }{Principal Components - if sizeshape=TRUE, the first variable of the PCs is size information (as log transformed Centroid Size)}
  \item{PCsym }{Principal Components of the Symmetrical Component}
  \item{PCasym }{Principal Components of the Asymmetrical Component}
  \item{PCscores }{PC scores}
  \item{PCscore_sym }{PC scores of the Symmetrical Component}
  \item{PCscore_asym }{PC scores of the Asymmetrical Component}
  \item{eigenvalues }{eigenvalues of the Covariance matrix}
  \item{eigensym }{eigenvalues of the "Symmetrical" Covariance matrix}
  \item{eigenasym }{eigenvalues of the "Asymmetrical" Covariance matrix}
  \item{Variance }{Table of the explained Variance by the PCs}
  \item{SymVar }{Table of the explained "Symmetrical" Variance by the PCs}
  \item{AsymVar }{Table of the explained "Asymmetrical" Variance by the PCs}
  \item{orpdata }{k x m x n array of the rotated configurations projected into tangent space}
  \item{rho }{vector of Riemannian distance from the mean}
}
\references{
  Dryden IL, and Mardia KV. 1998. Statistical shape analysis. Chichester.
  
Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape analysis of
symmetric structures: quantifying variation among individuals and
asymmetry. Evolution 56(10):1909-1920.

Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. «Semilandmarks in Three Dimensions,» in Modern Morphometrics in Physical Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer Academic/Plenum Publishers.


}
\author{
Stefan Schlager
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
  

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
