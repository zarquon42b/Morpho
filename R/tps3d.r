#' thin plate spline mapping
#' 
#' maps a datamatrix via thin plate spline between calculated by a reference on
#' a target configuration in 2D and 3D
#' 
#' 
#' @param M datamatrix - e.g. the matrix information of vertices of a given
#' surface
#' @param refmat reference matrix - e.g. landmark configuration on a surface
#' @param tarmat target matrix - e.g. landmark configuration on a target
#' surface
#' 
#' @param lambda integer: regularisation parameter of the TPS.
#' @return returns the warped datamatrix
#' @author Stefan Schlager
#' @seealso \code{\link{warp.mesh}}
#' @references Bookstein FL. 1989. Principal Warps: Thin-plate splines and the
#' decomposition of deformations. IEEE Transactions on pattern analysis and
#' machine intelligence 11(6).
#' @examples
#' 
#' require(Morpho)
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
#' }
#' 
#' @export
tps3d <- function(M,refmat,tarmat,lambda=1e-5)
{   
    q <- dim(M)[1]
    p <- dim(refmat)[1]
    m <- dim(refmat)[2]
    Lall <- CreateL(refmat,lambda=lambda, output="Linv")
    Linv <- Lall$Linv
    m2 <- rbind(tarmat,matrix(0,m+1,m))
    coeff <- matrix(NA,p+m+1,m)
    transM <- matrix(NA,q,m)
    coeff <- as.matrix(Linv%*%m2)
    #checks <- unlist(lapply(list(refmat,M,coeff), function(x){ out <- is.matrix(x) && is.numeric(x);return(out)}))
    #if (!prod(checks))
    #    stop("M, refmat and tarmat must be numeric matrices")
    transM <- .fx(refmat,M,coeff)
 
    return(transM)
    
}

  
