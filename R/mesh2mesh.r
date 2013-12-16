#' projects the vertices of a mesh onto the surface of another one.
#' 
#' projects the vertices of a mesh onto the surface of another one by searching
#' for the closest point (Euclidean distance or along vertex normals) on the
#' target by for each vertex.
#' 
#' needs trimesh_project and rayproject from trimesh-tools to be installed (
#' \url{http://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/})
#' 
#' @title projects the vertices of a mesh onto the surface of another one.
#' @param mesh1 mesh to project. Can be an object of class "mesh3d" or path to
#' an external mesh file (ply, obj, stl).
#' @param tarmesh mesh to project onto. Can be an object of class "mesh3d" or
#' path to an external mesh file (ply, obj, stl).
#' @param clean logical: request removing of dumpfiles.
#' @param cloud logical: if TRUE and mesh1 is an external mesh file, face
#' information will not be read from this file (saves time), as there is none
#' present.
#' @param sign logical: if TRUE, signed distances are returned.
#' @param tol numeric: maximum distance to search along ray, closest Euclidean
#' distance will be used, if tol is exceeded.
#' @param angmax numeric: maximum angle (in radians) of normals of original and
#' hit point are allowed to differ.
#' @param outname character: set name of dumpfile used in the process.
#' @param readback logical: whether the data is to be read into workspace.
#' @param inbound inverse search direction along rays.
#' @param strict logical: writes the value 1e12 into vertex quality if no face
#' is hit by the ray.
#' @param ignore.stdout suppress command line output.
#' @param mindist search in both directions of the ray and use closest points.
#' @return returns projected mesh
#' @author Stefan Schlager
#' @seealso \code{\link{ply2mesh}}, \code{\link{closemeshKD}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling.
#' @rdname mesh2mesh
#' @export
mesh2mesh <- function(mesh1,tarmesh,clean=TRUE,cloud=FALSE,sign=FALSE)
{
  options <- NULL
  if (sign)
    {options <- paste(options,"--sign")
   }
  it <- NULL
  if (!is.character(tarmesh))
    {mesh2ply(tarmesh,"target")
     tarply <- "target.ply"
   }
  else
    {tarply <- tarmesh
   }
  if (!is.character(mesh1))
    {
      vert <- t(mesh1$vb[1:3,])
      projBack(vert,tarply)
      it <- mesh1$it
    }
  
  else
    {
      if (!cloud)
        {
          it <- file2mesh(mesh1)$it
        }
      system(paste("trimesh_project ",mesh1," ",tarply,options, sep=""))
    }
  outmesh <- ply2mesh("out_cloud.ply",readnormals=TRUE)
  class(outmesh) <- c("mesh3d","shapes3d")
  outmesh$vb <- rbind(outmesh$vb,1) 
  outmesh$it <- it
  if (!is.null(it))
    {
      outmesh <- updateNormals(outmesh)
    }
  if (clean &&!is.character(tarmesh) )
    {
      if (!is.character(tarmesh))
        {unlink(c("target.ply"))
       }
      unlink("out_cloud.ply")
    }
  return(outmesh)
}
