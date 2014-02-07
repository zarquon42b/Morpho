\encoding{utf8}
\name{slider3d}
\alias{slider3d}
\title{slides Semilandmarks along curves and surfaces in 3D by minimising bending
energy of a thin-plate spline deformation.}
\usage{
slider3d(dat.array, SMvector, outlines = NULL, surp = NULL,
  sur.path = "sur", sur.name = NULL, meshlist = NULL, ignore = NULL,
  sur.type = "ply", tol = 1e-05, deselect = FALSE, inc.check = TRUE,
  recursive = TRUE, iterations = 0, initproc = TRUE, speed = TRUE,
  pairedLM = 0, weights = NULL, mc.cores = parallel::detectCores(),
  fixRepro = TRUE)
}
\arguments{
  \item{dat.array}{Input k x m x n real array, where k is
  the number of points, m is the number of dimensions, and
  n is the sample size. Ideally the dimnames[[3]] vector
  contains the names of the surface model (without file
  extension) - e.g. if the model is named "surface.ply",
  the name of the corresponding matrix of the array would
  be "surface"}

  \item{SMvector}{A vector containing the landmarks on the
  curve(s) and surfaces that are allowed to slide}

  \item{outlines}{A vector (or if threre are several
  curves) a list of vectors (containing the rowindices) of
  the (Semi-)landmarks forming the curve(s) in the
  successive position on the curve - including the
  beginning and end points, that are not allowed to slide.}

  \item{surp}{A vector containing Semilandmarks positioned
  on surfaces.}

  \item{sur.path}{Path to the surface models (e.g. ply,
  obj, stl files)}

  \item{sur.name}{character vector: containing the
  filenames of the corresponding surfaces - e.g. if the
  dat.array[,,i] belongs to surface_i.ply, sur.name[i]
  would be surface_i.ply. Only necessary if dat.array does
  not contain surface names.}

  \item{meshlist}{list containing triangular meshes of
  class 'mesh3d', for example imported with
  \code{\link{mesh2ply}} or \code{\link{file2mesh}} in the
  same order as the specimen in the array (see examples
  below)}

  \item{ignore}{vector containing indices of landmarks that
  are to be ignored. Indices of outlines/surfaces etc will
  be updated automatically.}

  \item{sur.type}{character:if all surfaces are of the same
  file format and the names stored in dat.array, the file
  format will be specified here.}

  \item{tol}{numeric: Threshold for convergence in the
  sliding process}

  \item{deselect}{Logical: if TRUE, the SMvector is
  interpreted as those landmarks, that are not allowed to
  slide.}

  \item{inc.check}{Logical: if TRUE, the program stops when
  convergence criterion starts increasing and reports
  result from last iteration.}

  \item{speed}{Logical: if TRUE, only a partial procrustes
  fit will be performed - this is faster and can be
  required, when large samples are processed.}

  \item{recursive}{Logical: if TRUE, during the iterations
  of the sliding process, the outcome of the previous
  iteration will be used.  Otherwise the original
  configuration will be used in all iterations.}

  \item{iterations}{integer: select manually the max.
  number of iterations that will be performed during the
  sliding process (usefull, when there is very slow
  convergence).  0 means iteration until convergence.}

  \item{initproc}{requests initial Procrustes fit before
  sliding.}

  \item{pairedLM}{A X x 2 numeric matrix with the indices
  of the rows containing paired Landmarks. E.g. the left
  column contains the lefthand landmarks, while the right
  side contains the corresponding right hand landmarks. -
  This will ideally create symmetric mean to get rid of
  assymetry.}

  \item{weights}{vector: assign a weight to each landmark:
  the smaller the value is, the less it will be affected by
  sliding. 0 = fix. This is highly experimental!!!}

  \item{mc.cores}{integer: determines how many cores to use
  for the computation. The default is autodetect. But in
  case, it doesn't work as expected cores can be set
  manually. In Windows, parallel processing is disabled.}

  \item{fixRepro}{logical: if \code{TRUE}, fix landmarks
  will also be projected onto the surface. If you have
  landmarks not on the surface, select
  \code{fixRepro=FALSE}}
}
\value{
\item{dataslide }{array containing slidden Landmarks in the
original space - not yet processed by a Procrustes
analysis} \item{vn.array }{array containing landmark
normals}
}
\description{
slides Semilandmarks along curves and surfaces in 3D. The
positions on the surface are sought which minimise bending
energy (of a thin-plate spline deformation)
}
\section{Warning}{
  Depending on the size of the suface meshes and especially
  the amount of landmarks this can use an extensive amount
  of your PC's resources, especially when running in
  parallel. As the computation time and RAM usage of matrix
  algebra involved is quadratic to the amount of landmarks
  used, doubling the amount of semi-landmarks will
  quadruple computation time and system resource usage. You
  can easily stall you computer with this function with
  inappropriate data.
}
\examples{
\dontrun{
data(nose)
###create mesh for longnose
longnose.mesh <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
### write meshes to disk
mesh2ply(shortnose.mesh, filename="shortnose")
mesh2ply(longnose.mesh, filename="longnose")

## create landmark array
data <- bindArr(shortnose.lm, longnose.lm, along=3)
dimnames(data)[[3]] <- c("shortnose", "longnose")

# define fix landmarks
fix <- c(1:5,20:21)
# define surface patch by specifying row indices of matrices
# all except those defined as fix
surp <- c(1:nrow(shortnose.lm))[-fix]

slide <- slider3d(data, SMvector=fix, deselect=TRUE, surp=surp,
                  sur.path=".",iterations=1)
                  # sur.path="." is the current working directory

# now one example with meshes in workspace
## to reduce this example's computation time,
# we only use the first 50 right hand semi-landmarks
surp <- surp[1:50]
meshlist <- meshlist <- list(shortnose.mesh,longnose.mesh)

slide <- slider3d(data[1:57,,], SMvector=fix, deselect=TRUE, surp=surp,
                  sur.path=".",iterations=1, meshlist=meshlist,
                  mc.cores=1,fixRepro=FALSE)
require(rgl)
## visualize sliding
deformGrid3d(slide$dataslide[,,1],shortnose.lm,ngrid = 0)
## these are fix
spheres3d(slide$dataslide[fix,,1],col=4,radius=0.7)
}
}
\author{
Stefan Schlager
}
\references{
Klingenberg CP, Barluenga M, and Meyer A. 2002. Shape
analysis of symmetric structures: quantifying variation
among individuals and asymmetry. Evolution
56(10):1909-1920.

Gunz, P., P. Mitteroecker, and F. L. Bookstein. 2005.
Semilandmarks in Three Dimensions, in Modern Morphometrics
in Physical Anthropology. Edited by D. E. Slice, pp. 73-98.
New York: Kluwer Academic/Plenum Publishers.

Schlager S. 2012. Sliding semi-landmarks on symmetric
structures in three dimensions. American Journal of
Physical Anthropology, 147(S52):261. URL:
http://dx.doi.org/10.1002/ajpa.21502.

Schlager S. 2013. Soft-tissue reconstruction of the human
nose: population differences and sexual dimorphism. PhD
thesis,
\enc{Universitätsbibliothek}{Universitaetsbibliothek}
Freiburg.  URL:
\url{http://www.freidok.uni-freiburg.de/volltexte/9181/}.
}
\seealso{
\code{\link{relaxLM}}
}

