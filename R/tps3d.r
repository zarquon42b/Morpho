#' thin plate spline mapping (2D and 3D) for coordinates and triangular meshes
#' 
#' maps landmarks or a triangular mesh via thin plate spline based on a reference and
#' a target configuration in 2D and 3D
#' 
#' 
#' @param x matrix - e.g. the matrix information of vertices of a given
#' surface or a triangular mesh of class "mesh3d"
#' @param refmat reference matrix - e.g. landmark configuration on a surface
#' @param tarmat target matrix - e.g. landmark configuration on a target
#' surface
#' 
#' @param lambda numeric: regularisation parameter of the TPS.
#' @param ... additional arguments, currently not used.
#' @return returns the deformed input
#' @author Stefan Schlager
#' @seealso \code{\link{computeTransform}, \link{applyTransform}}
#' @references Bookstein FL. 1989. Principal Warps: Thin-plate splines and the
#' decomposition of deformations. IEEE Transactions on pattern analysis and
#' machine intelligence 11(6).
#' @examples
#' 
#' data(nose)
#' ## define some landmarks
#' refind <- c(1:3,4,19:20)
#' ## use a subset of shortnose.lm as anchor points for a TPS-deformation
#' reflm <- shortnose.lm[refind,]
#' tarlm <- reflm
#' ##replace the landmark at the tip of the nose with that of longnose.lm
#' tarlm[4,] <- longnose.lm[4,]
#' ##  deform a set of semilandmarks by applying a TPS-deformation
#' ##  based on 5 reference points
#' deformed <- tps3d(shortnose.lm, reflm, tarlm)
#' \dontrun{
#' ##visualize results by applying a deformation grid
#' deformGrid3d(shortnose.lm,deformed,ngrid = 5)
#'
#' 
#' data(nose)##load data
#' ##warp a mesh onto another landmark configuration:
#' warpnose.long <- tps3d(shortnose.mesh,shortnose.lm,longnose.lm)
#' 
#' 
#' require(rgl)
#' shade3d(warpnose.long,col=skin1)
#' }
#' 
#' data(boneData)
#' ## deform mesh belonging to the first specimen
#' ## onto the landmark configuration of the 10th specimen
#'
#' \dontrun{
#' warpskull <- tps3d(skull_0144_ch_fe.mesh,boneLM[,,1],
#'                      boneLM[,,10])
#' ## render deformed mesh and landmarks
#' shade3d(warpskull, col=2, specular=1)
#' spheres3d(boneLM[,,1])
#' ## render original mesh
#' shade3d(skull_0144_ch_fe.mesh, col=3, specular=1)
#' spheres3d(boneLM[,,10])
#' 
#' }
#' @export
tps3d <- function(x,refmat,tarmat,lambda=1e-8,...) {
    coeff <- computeTransform(x=tarmat,y=refmat,lambda=lambda,type="tps")
    transM <- applyTransform(x,coeff)
    return(transM)
    
}
