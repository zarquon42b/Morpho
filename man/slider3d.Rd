\name{slider3d}

\alias{slider3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 slides Semilandmarks along curves and surfaces in 3D by minimising
 bending energy of a thin-plate spline deformation.

}
\description{slides Semilandmarks along curves and surfaces in 3D. The
  positions on the surface are sought which minimise bending energy (of
  a thin-plate spline deformation)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
slider3d(dat.array, SMvector, outlines=NULL, surp=NULL, sur.path="sur",
sur.name=NULL, ignore=NULL, sur.type="ply", tol=1e-05, deselect=FALSE,
inc.check=TRUE, recursive=TRUE, iterations=0, initproc=TRUE, speed=TRUE,
pairedLM=0, weights=NULL, mc.cores=detectCores())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat.array}{Input k x m x n real array, where k is the number of points, m is the number of dimensions, and n is the sample size. Ideally the dimnames[[3]] vector contains the names of the surface model (without file extension) - e.g. if the model is named "surface.ply", the name of the corresponding matrix of the array would be "surface"
%%     ~~Describe \code{dat.array} here~~
}
  \item{SMvector}{A vector containing the landmarks on the curve(s) and
    surfaces that are allowed to slide
%%     ~~Describe \code{SMvector} here~~
}
  \item{outlines}{A vector (or if threre are several curves) a list of vectors (containing the rowindices) of the (Semi-)landmarks forming the curve(s) in the successive position on the curve - including the beginning and end points, that are not allowed to slide.
%%     ~~Describe \code{outlines} here~~
}
  \item{surp}{A vector containing Semilandmarks positioned on surfaces.
%%     ~~Describe \code{surp} here~~
}
  \item{sur.path}{Path to the surface models (e.g. ply, obj, stl files)
%%     ~~Describe \code{sur.path} here~~
}
  \item{sur.name}{character vector: containing the filenames of the corresponding surfaces - e.g. if the dat.array[,,i] belongs to surface_i.ply, sur.name[i] would be surface_i.ply. Only necessary if dat.array does not contain surface names. 
%%     ~~Describe \code{sur.name} here~~
}
 \item{ignore}{vector containing rownames of landmarks that are to be ignored. Indices of outlines/surfaces etc will be updated automatically
}
  \item{sur.type}{character:if all surfaces are of the same file format and the names stored in dat.array, the file format will be specified here.
%%     ~~Describe \code{sur.type} here~~
}
  \item{tol}{numeric: Threshold for convergence in the sliding process
%%     ~~Describe \code{tol} here~~
}
  \item{deselect}{Logical: if TRUE, the SMvector is interpreted as those landmarks, that are not allowed to slide.
%%     ~~Describe \code{deselect} here~~
}
  \item{inc.check}{Logical: if TRUE, the program stops when convergence criterion starts increasing and reports result from last iteration.
%%     ~~Describe \code{deselect} here~~
}
 \item{speed}{Logical: if TRUE, only a partial procrustes fit will be performed  - this is faster and can be required, when large samples are processed.
%%     ~~Describe \code{deselect} here~~
}
  \item{recursive}{Logical: if TRUE, during the iterations of the sliding process, the outcome of the previous iteration will be used. 
                    Otherwise the original configuration will be used in all iterations. 
%%     ~~Describe \code{recursive} here~~
}
  \item{iterations}{integer: select manually the max. number of iterations that will be performed during the sliding process 
                    (usefull, when there is very slow convergence). 
                              0 means iteration until convergence.
%%     ~~Describe \code{iterations} here~~
}
 
  \item{initproc}{requests initial Procrustes fit before sliding.
%%     ~~Describe \code{initproc} here~~
}
  \item{pairedLM}{A X x 2 numeric matrix with the indices of the rows containing paired Landmarks. E.g. the left column contains the lefthand landmarks, while the right side contains the corresponding right hand landmarks. - This will ideally create symmetric mean to get rid of assymetry along outline after first relaxation. The effect is not yet tested.
%%     ~~Describe \code{pairedLM} here~~
}  
\item{weights}{vector: assign a weight to each landmark: the smaller the
  value is, the less it will be affected by sliding. 0 = fix. This is
  highly experimental!!!
%%     ~~Describe \code{initproc} here~~
}
\item{mc.cores}{integer: determines how many cores to use for the
  computation. The default is autodetect. But in case, it doesn't work
  as expected cores can be set manually. In Windows, parallel processing
  is disabled.
%%     ~~Describe \code{initproc} here~~
}
}
\details{needs trimesh_project to be installed (\url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/})

%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{dataslide }{array containing slidden Landmarks in the original
   space - not yet processed by a Procrustes analysis.}
 \item{vn.array }{array containing Landmark normals}
%% ...
}
\references{Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape
  analysis of symmetric structures: quantifying variation among
  individuals and asymmetry. Evolution 56(10):1909-1920.
  

  Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005. Semilandmarks in Three Dimensions, in Modern Morphometrics in Physical Anthropology. Edited by D. E. Slice, pp. 73-98. New York: Kluwer Academic/Plenum Publishers.
%% ~put references to the literature/web site here ~
}
\author{Stefan Schlager
%%  ~~who you are~~
}


\seealso{
\code{\link{relaxLM}}
}

