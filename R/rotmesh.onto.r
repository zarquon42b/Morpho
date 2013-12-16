#' rotate ,scale and translate a mesh based on landmark information.
#' 
#' rotates and reflects a mesh onto by calculating the transformation from two
#' sets of referenced landmarks.
#' 
#' 
#' @param mesh object of class mesh3d.
#' @param refmat k x m matrix with landmarks on the mesh
#' @param tarmat k x m matrix as target configuration
#' @param adnormals logical - if TRUE, vertex normals will be recomputed after
#' rotation. If \code{mesh} has normals and adnormals=FALSE, the existing
#' normals are rotated by the same rotation matrix as the mesh's vertices.
#' @param scale logical: if TRUE the mesh will be scaled according to the size
#' of the target.
#' @param reflection logical: allow reflection.
#' 
#' @return
#' \item{mesh }{rotated mesh}
#' \item{yrot }{rotated refmat}
#' @author Stefan Schlager
#' @seealso \code{\link{file2mesh}},\code{\link{warp.mesh}}
#' ,\code{\link{rotonto}},\code{\link{mesh2ply}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(boneData)
#' ## rotate, translate and scale the mesh belonging to the first specimen
#' ## onto the landmark configuration of the 10th specimen
#' rotmesh <- rotmesh.onto(skull_0144_ch_fe.mesh,boneLM[,,1],
#'                         boneLM[,,10], scale=TRUE)
#' \dontrun{
#' ## render rotated mesh and landmarks
#' shade3d(rotmesh$mesh, col=2, specular=1)
#' spheres3d(boneLM[,,1])
#' ## render original mesh
#' shade3d(skull_0144_ch_fe.mesh, col=3, specular=1)
#' spheres3d(boneLM[,,10])
#' }
#' 
#' @export
rotmesh.onto <- function(mesh, refmat, tarmat, adnormals=FALSE, scale=FALSE, reflection=FALSE)
{
  rot <- rotonto(tarmat,refmat,scale=scale,reflection=reflection)
  mesh$vb[1:3,] <- mesh$vb[1:3,]-rot$transy
  mesh$vb[1:3,] <- t(t(mesh$vb[1:3,])%*%rot$gamm)
  if (scale)
      mesh$vb[1:3,] <- mesh$vb[1:3,]*rot$bet
    
  mesh$vb[1:3,] <- mesh$vb[1:3,]+rot$trans
  if (sign(det(rot$gamm) < 0 && reflection))
      mesh <- conv2backf(mesh)
  if (adnormals) 
      mesh <- updateNormals(mesh)
  if (!is.null(mesh$normals) && !adnormals)
      mesh$normals[1:3,] <- t(t(mesh$normals[1:3,]) %*% rot$gamm)
   
  return(list(mesh=mesh,yrot=rot$yrot))
}
