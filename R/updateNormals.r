#' Compute face or vertex normals of a triangular mesh
#' 
#' Compute face or vertex normals of a triangular mesh of class "mesh3d"
#' 
#' 
#' @param x triangular mesh of class "mesh3d"
#' @param angle logical: if TRUE, angle weighted normals are used.
#' @return \code{updateNormals} returns mesh with updated vertex normals.
#' 
#' \code{facenormals} returns an object of class "mesh3d" with
#' \item{vb }{faces' barycenters}
#' \item{normals }{faces' normals}
#' @note only supports triangular meshes
#' @author Stefan Schlager
#' @seealso \code{\link{ply2mesh}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling, .
#' 
#' @examples
#' 
#' require(rgl)
#' require(Morpho)
#' data(nose)
#' ### calculate vertex normals
#' shortnose.mesh$normals <- NULL ##remove normals
#' \dontrun{
#' shade3d(shortnose.mesh,col=3)##render
#' }
#' shortnose.mesh <- updateNormals(shortnose.mesh)
#' \dontrun{
#' rgl.clear()
#' shade3d(shortnose.mesh,col=3)##smoothly rendered now
#' }
#' ## calculate facenormals
#' facemesh <- facenormals(shortnose.mesh)
#' \dontrun{
#' plotNormals(facemesh,long=0.01)
#' points3d(vert2points(facemesh),col=2)
#' wire3d(shortnose.mesh)
#' }
#' @rdname updateNormals
#' @export
updateNormals <- function(x,angle=TRUE) 
{
    if (!is.null(x$ib)) {
        x <- quad2trimesh(x,FALSE)
        message("mesh was converted to triangular mesh")
    }
    vb <- x$vb
    ## Make sure v is homogeneous with unit w
    if (nrow(vb) == 3)
        vb <- rbind(vb,1)
    else
        vb <- t(t(vb)/vb[4,])
    vb <- vb[1:3,]
    if (!is.matrix(vb) || !is.numeric(vb))
        stop("vertices must be a numeric matrix")
    if (!is.null(x$it)) {
        if (!is.matrix(x$it) || !is.numeric(x$it))
            stop("faces indices must be stored as numeric matrix")
        rangeit <- range(x$it)
        if (rangeit[1] < 1 || rangeit[2] > ncol(vb))
            stop("faces point beyond range of vertices")
        it <- x$it-1
        out <- .Call("updateVertexNormals",vb,it,angle)
        normals <- rbind(out,1)
        x$normals <- normals
    } else {
        message("mesh has no triangular faces, nothing to be done")
    }
    return(x)
}
#' @rdname updateNormals
#' @export
facenormals <- function(x) 
{
    barymesh <- list()
    barymesh$vb <- rbind(t(barycenter(x)),1)    
    vb <- x$vb
    ## Make sure v is homogeneous with unit w
    if (nrow(vb) == 3)
        vb <- rbind(vb, 1)
    else
        vb <- t( t(vb)/vb[4,] )
    vb <- vb[1:3,]
    if (!is.matrix(vb) || !is.numeric(vb))
        stop("vertices must be a numeric matrix")
    if (!is.null(x$it)) {
        it <- x$it-1
        
    }
    else
        stop("mesh has no triangular faces")

    out <- .Call("updateFaceNormals",vb,it)
    normals <- out
    class(barymesh) <- "mesh3d"
    barymesh$normals <- normals
    
    return(barymesh)
}


