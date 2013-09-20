\name{placePatch}
\alias{placePatch}
\title{
 Project semi-landmarks from a predefined atlas onto all specimen in a sample
}
\description{
 Project semi-landmarks from a predefined atlas onto all specimen in a
 sample. Various mechanisms are implemented to avoid errorneous
 placement on the wrong surface layer (e.g. inside the bone). 
}
\usage{
placePatch(atlas, dat.array, path, prefix = NULL, fileext = ".ply",
           ray = TRUE, inflate = NULL, tol=inflate, relax.patch = TRUE,
           keep.fix = NULL, rhotol = NULL)
}
\arguments{
  \item{atlas}{object of class "atlas" created by \code{\link{createAtlas}}
  }
  \item{dat.array}{k x 3 x n array containing reference landmarks of the
    sample or a matrix in case of only one target specimen.
  }
  \item{path}{character: specify the directory where the surface meshes
    of the sample are stored. In case \code{dat.array} is a matrix,
    \code{path} is needed to specifiy not only the directory, but the
   specific surface mesh file (without fileextension if specified by
   \code{fileext}). See examples below.
  }
  \item{prefix}{character: append a prefix to the specimens names
    (stored in \code{dimnames(dat.array)[[3]]}) to match the
    corresponding file names.
  }
  \item{fileext}{character: file extension of the surface meshes.
  }
  
  \item{ray}{logical: projection will be along surface normals instead
    of simple closest point search.
  }
  \item{inflate}{inflate (or deflate - if negative sign) the
    semilandmarks along the normals of the deformed atlas to make sure
    that they stay on the outside (inside) of the target mesh.
  }
  \item{tol}{numeric: threshold to follow the ray back after
    inflation. See details below. If no surface is hit after \code{tol}
    mm, the simple closest point will be used.
  }
  \item{relax.patch}{logical: request relaxation minimising bending
    energy toward the atlas.
  }
  \item{keep.fix}{integer: rowindices of those landmarks that are not
    allowed to be relaxed in case \code{relax.patch=TRUE}. If not
    specified, all landmarks will be kept fix.
  }
  \item{rhotol}{numeric: maximum amount of deviation a hit point's
    normal is allowed to deviate from the normal defined on the
    atlas. If \code{relax.patch=TRUE}, those points exceeding this value
    will be relaxed freely (i.e. not restricted to tangent plane).
  }
}
\details{This function allows the (relatively) easy projection of
  surface points defined on an atlas onto all surface of a given sample
  by Thin-Plate Spline deformation and additional mechanisms to avoid
  distortions. The algorithm can be outlined as followed.
  \enumerate{
  \item{deform atlas onto targets by TPS based on predefined
    landmarks (and curves).
  }
  \item{project coordinates on deformed atlas onto target mesh
  }
  \item{'inflate' or 'deflate' configuration along their normals to make
    sure all coordinates are on the outside/inside 
  }
  \item{Project inflated points back onto surface along these
    normals.
  }
  \item{Check if normals are roughly pointing into the same direction
    as those on the (deformed) atlas.
  }
   \item{Relax all points against atlas.
   }
 }
 
}
\value{array containing the projected coordinates appended to the
  data.array specified in the input. In case dat.array is a matrix only
  a matrix is returned.
}
\references{mydiss

  
}
\author{
Stefan Schlager}
\note{needs additional command line tools "trimesh-tools" installed (\url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/}).

}

\examples{
data(nose)
require(rgl)
###create mesh for longnose
longnose.mesh <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
## create atlas
fix <- c(1:5,20:21)
atlas <- createAtlas(shortnose.mesh, landmarks =
           shortnose.lm[fix,], patch=shortnose.lm[-c(1:5,20:21),])
## view atlas
plotAtlas(atlas)

## create landmark array with only fix landmarks
data <- bindArr(shortnose.lm[fix,], longnose.lm[fix,], along=3)
dimnames(data)[[3]] <- c("shortnose", "longnose")

### write meshes to disk
mesh2ply(shortnose.mesh, filename="shortnose")
mesh2ply(longnose.mesh, filename="longnose")
\dontrun{
patched <- placePatch(atlas, data, path="./", inflate=5)
## now browse through placed patches
checkLM(patched, path="./", atlas=atlas)
}
## same example with only one target specimen
data <- longnose.lm[fix, ]
\dontrun{
patched <- placePatch(atlas, data, path="longnose", inflate=5)
wire3d(longnose.mesh,col=3)
spheres3d(patched)
}

}
\seealso{
 \code{\link{createAtlas}, \link{relaxLM}, \link{checkLM}, \link{slider3d}, \link{warp.mesh}}
}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
