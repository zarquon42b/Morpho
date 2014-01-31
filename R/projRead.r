projBack <- function(data,surface,dataname=NULL,outname=NULL,smooth=TRUE,ignore.stdout=FALSE,sign=FALSE)
{
  checkCLI <- try(system("trimesh_project",intern=TRUE),silent=TRUE)
   if(class(checkCLI) == "try-error")
      stop("please install trimesh-tools")
  options <- NULL
  
  if (!smooth)
      options <- paste(options,"--nosmooth")
  if (sign)
      options <- paste(options,"--sign")
  if (is.null(dataname))
      dataname <- "out"
  write.obj(cbind("v",data),filename=dataname)
  if (is.null(outname))
      command <- paste("trimesh_project"," ",dataname,".obj"," ",surface,options,sep="")  
  else
      command <- paste("trimesh_project"," ",dataname,".obj"," ",surface," -o ",outname,options,sep="")
  
  system(command,ignore.stdout=ignore.stdout)
  unlink(paste(dataname,".obj",sep="")) #clean up
}



#' Project points onto the closest point on a mesh
#' 
#' project points onto a given surface and return projected points and normals.
#' 
#' 
#' @param lm m x 3 matrix containing 3D coordinates.
#' @param mesh character: specify path to mesh file.
#' @param readnormals logical: return normals of projected points.
#' @param smooth logical: rerturn smoothed normals.
#' @param sign logical: request signed distances.
#' @param \dots additional arguments not used at the moment.
#' @return if readnormals = FALSE, a m x 3 matrix containing projected points
#' is returned, otherwise a list, where
#' \item{vb }{3 x m matrix containing projected points}
#' \item{normals }{3 x m matrix containing normals}
#' \item{quality }{vector containing distances }
#' @note The usage of this function requires the command line tools from
#' trimesh-tools
#' (https://sourceforge.net/projects/morpho-rpackage/files/Auxiliaries/)
#' installed.
#' @author Stefan Schlager
#' @seealso \code{\link{closemeshKD}}
#' @references Detection of inside/outside uses the algorithm proposed in:
#' 
#' Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance
#' Fields From Triangle Meshes. Informatics and Mathematical Modelling.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' data(nose)
#' \dontrun{
#' repro <- projRead(shortnose.lm,shortnose.mesh)
#' }
#' 
#' @importFrom Rvcg vcgClost vcgImport
#' @export
projRead <- function(lm, mesh,readnormals=TRUE, smooth=TRUE, sign=TRUE,...)
{
    if (is.character(mesh))
        mesh <- vcgImport(mesh)
    
    data <- vcgClost(lm, mesh, smoothNormals=smooth,sign=sign,borderchk=FALSE)
    if (!readnormals)
        data <- vert2points(data)
    return(data)
}
