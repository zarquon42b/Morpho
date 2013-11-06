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
#' @param clean logical: remove dumpfiles.
#' @param smooth logical: rerturn smoothed normals.
#' @param ignore.stdout logical: supress console messages from system calls.
#' @param sign logical: request signed distances.
#' @param prodump character: name of the dumpfile storing the projected points'
#' coordinates (useful when using parallel backend).
#' @param lmdump character: name of the dumpfile storing data to be projected
#' (useful when using parallel backend).
#' @return if readnormals = FALSE, a m x 3 matrix containing projected points
#' is returned, otherwise a list, where
#' \item{vb }{3 x m matrix containing projected points}
#' \item{normals }{3 x m matrix containing normals}
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
#' 
#' @export projRead
projRead <- function(lm,mesh,readnormals=TRUE,clean=TRUE,smooth=TRUE,ignore.stdout=FALSE,sign=FALSE,lmdump=NULL,prodump=NULL)
{

  checkCLI <- try(system("trimesh_project",intern=TRUE),silent=TRUE)
   if(class(checkCLI) == "try-error")
      stop("please install trimesh-tools")
  if (is.null(prodump))
    prodump <- "out_cloud.ply"
  if (is.null(lmdump))
    lmdump <- "out"
  if (is.character(mesh))
      projBack(lm,mesh,ignore.stdout=ignore.stdout,sign=sign,smooth=smooth,outname=prodump,dataname=lmdump)
  else {
      dumpname <- paste(prodump,"dump0",sep="")
      dumpfile <- paste(dumpname,".ply",sep="")
      mesh2ply(mesh,dumpname)
      projBack(lm,dumpfile,smooth=smooth,ignore.stdout=ignore.stdout,sign=sign,outname=prodump,dataname=lmdump)
      unlink(dumpfile)
  }
  
  data <- ply2mesh(prodump,readnormals=readnormals, silent=ignore.stdout)
  if (clean)
      unlink(prodump)
  
  return(data)
}
