#' warping a mesh onto another configuration
#' 
#' warps an the surface of a mesh3d object onto another configuration via
#' reference and target landmark configuration by using a thin-plate spline
#' interpolation.
#' 
#' the surface is mapped via the tps3d function onto the target shape.
#' 
#' @param mesh object of class "mesh3d"
#' 
#' @param matr matrix of landmarks on the reference surface
#' 
#' @param matt matrix of corresponding landmarks on the target surface
#' 
#' @param updateNormals Logical: requests the (re)calculation of vertex
#' normals.
#' @param lambda integer: regularisation parameter of the TPS.
#' @param silent logical: suppress messages.
#' @return
#' object of class "mesh3d"
#' @author Stefan Schlager
#' @seealso
#' \code{\link{ply2mesh},\link{file2mesh},\link{mesh2ply},\link{warpmovie3d},
#' \link{rotmesh.onto}}
#' @examples
#' 
#' require(rgl)
#' data(nose)##load data
#' ##warp a mesh onto another landmark configuration:
#' warpnose.long <- warp.mesh(shortnose.mesh,shortnose.lm,longnose.lm)
#'\dontrun{
#' shade3d(warpnose.long,col=skin1)
#' }
#' 
#' data(boneData)
#' ## deform mesh belonging to the first specimen
#' ## onto the landmark configuration of the 10th specimen
#' \dontrun{
#' warpskull <- warp.mesh(skull_0144_ch_fe.mesh,boneLM[,,1],
#'                      boneLM[,,10])
#' ## render deformed mesh and landmarks
#' shade3d(warpskull, col=2, specular=1)
#' spheres3d(boneLM[,,1])
#' ## render original mesh
#' shade3d(skull_0144_ch_fe.mesh, col=3, specular=1)
#' spheres3d(boneLM[,,10])
#' }
#' 
#' 
#' 
#' @export
warp.mesh <- function(mesh,matr,matt,lambda=0,updateNormals=TRUE, silent=FALSE)
{
    vert <- t(mesh$vb[1:3,])
    if (!silent)
        cat("calculating spline...\n")
    warp <- tps3d(vert,matr,matt,lambda=lambda)
    mesh$vb <- rbind(t(warp),1)
    mesh$normals <- NULL
    testref <- rotonto(matr,matt)$reflect
    if(testref == 1)
        mesh <- conv2backf(mesh)
    
    if(updateNormals) {
        if (!silent)
            cat("updating normals...\n")
        mesh <- updateNormals(mesh)
    }
    return(mesh)
}

