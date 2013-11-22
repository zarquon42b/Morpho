#' A toolbox providing methods for data-acquisitiopn, visualisation and
#' statistical methods related to Geometric Morphometrics and shape analysis
#' 
#' A toolbox for Morphometric calculations. Including sliding operations for
#' Semilandmarks, importing, exporting and manipulating of 3D-surface meshes
#' and semi-automated placement of surface landmarks. For gaining full
#' functionality, the command line tools provided on
#' \url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/}
#' have to be installed.
#' 
#' \tabular{ll}{
#' Package: \tab Morpho\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0-2.131122\cr
#' Date: \tab 2013-11-22\cr
#' License: \tab GPL\cr
#' LazyLoad: \tab yes\cr
#' }
#' 
#' @name Morpho-package
#' @aliases Morpho-package Morpho
#' @docType package
#' @note For full functionality please install trimesh-tools as described here:
#' \url{http://sourceforge.net/p/morpho-rpackage/wiki/Installation_Morpho/#installation-of-the-command-line-tools-strongly-recommended}
#'
#' The pdf-version of Morpho-help can be obtained from CRAN on \url{http://cran.r-project.org/web/packages/Morpho/Morpho.pdf}
#' @author Stefan Schlager
#' 
#' Maintainer: Stefan Schlager <stefan.schlager@@uniklinik-freiburg.de>
#' @references Schlager S. 2013. Soft-tissue reconstruction of the human nose:
#' population differences and sexual dimorphism. PhD thesis,
#' \enc{Universitätsbibliothek}{Universitaetsbibliothek} Freiburg. URL:
#' \url{http://www.freidok.uni-freiburg.de/volltexte/9181/}.
#' @encoding utf8
#' @keywords package
#' @useDynLib Morpho
#' @import doParallel
#' @importFrom colorRamps blue2green2red
#' @importFrom foreach foreach '%dopar%' '%do%'
#' @importFrom MASS ginv
#' @importFrom Matrix sparseMatrix diag crossprod solve
#' @importFrom parallel mclapply detectCores
#' @importFrom rgl lines3d open3d points3d rgl.bg rgl.bringtotop rgl.clear rgl.close  rgl.cur rgl.pop rgl.snapshot shade3d spheres3d text3d translate3d wire3d
#' @importFrom yaImpute ann
#' @importClassesFrom Matrix dgCMatrix dgeMatrix dsCMatrix dtCMatrix sparseMatrix
NULL


#' Landmarks and a triangular mesh
#' 
#' Landmarks on the osseous human nose and a triangular mesh representing this
#' structure.
#' 
#' 
#' @name boneData
#' @aliases boneLM skull_0144_ch_fe.mesh
#' @docType data
#' @format \code{boneLM}: A 10x3x80 array containing 80 sets of 3D-landmarks
#' placed on the human osseous nose.
#' 
#' \code{skull_0144_ch_fe.mesh}: The mesh representing the area of the first
#' individual of \code{boneLM}
#' @keywords datasets
NULL





#' predefined colors for bone and skin
#' 
#' predefined colors for bone and skin
#' 
#' available colors are:
#' 
#' bone1
#' 
#' bone2
#' 
#' bone3
#' 
#' skin1
#' 
#' skin2
#' 
#' skin3
#' 
#' skin4
#' 
#' @name colors
#' @export bone1 bone2 bone3 skin1 skin2 skin3 skin4
#' @aliases bone1 bone2 bone3 skin1 skin2 skin3 skin4
#' @docType data
#' @keywords datasets
NULL


#' landmarks and a triangular mesh representing a human nose
#' 
#' triangular mesh representing a human nose and two matrices containing
#' landmark data
#' 
#' 
#' @name nose
#' @aliases shortnose.mesh shortnose.lm longnose.lm
#' @docType data
#' @format \code{shortnose.mesh}: A triangular mesh of class 'mesh3d'.
#' 
#' \code{shortnose.lm}: matrix containing example landmark data placed on
#' \code{shortnose.mesh}.
#' 
#' \code{longnose.lm}: matrix containing example landmark data representing a
#' caricaturesquely deformed human nose.
#' @keywords datasets
NULL

## document deprecated functions
#'  @title deprecated functions of Morpho
#' @name deprecated
#' @rdname Morpho-deprecated
#' @keywords internal
NULL


#' @rdname Morpho-deprecated
#' @export deform.grid
deform.grid <- function (...)
{
  .Deprecated("deformGrid3d", package="Morpho")
  deformGrid3d(...)
}

#' @rdname Morpho-deprecated
#' @export regdist.raw
regdist.raw <- function (...)
{
  .Deprecated("regdist", package="Morpho")
  regdist(...)
}
