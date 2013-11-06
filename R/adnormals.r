#' Compute face or vertex normals of a triangular mesh
#' 
#' Compute face or vertex normals of a triangular mesh of class "mesh3d"
#' 
#' 
#' @param x triangular mesh of class "mesh3d"
#' @param angle logical: if TRUE, angle weighted normals are used.
#' @return adnormals returns mesh with updated vertex normals.
#' 
#' facenormals returns an object of class "mesh3d". With
#' \item{vb }{faces' barycenters}
#' \item{normals }{faces' normals}
#' @note only supports triangular meshes
#' @author Stefan Schlager
#' @seealso \code{\link{ply2mesh}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling, .
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' ### calculate vertex normals
#' shortnose.mesh$normals <- NULL ##remove normals
#' \dontrun{
#' shade3d(shortnose.mesh,col=3)##render
#' }
#' shortnose.mesh <- adnormals(shortnose.mesh)
#' \dontrun{
#' rgl.clear()
#' shade3d(shortnose.mesh,col=3)##smoothly rendered now
#' }
#' ## calculate facenormals
#' facemesh <- facenormals(shortnose.mesh)
#' plotNormals(facemesh,long=0.01)
#' \dontrun{
#' points3d(vert2points(facemesh),col=2)
#' wire3d(shortnose.mesh)
#' }
#' @rdname adnormals
#' @export adnormals
adnormals <- function(x,angle=TRUE) 
{
    v <- x$vb
    ## Make sure v is homogeneous with unit w
    if (nrow(v) == 3)
        v <- rbind(v, 1)
    else
        v <- t( t(v)/v[4,] )
    normals <- v*0
    v <- v[1:3,]
    it <- x$it
    storage.mode(v) <- "double"
    storage.mode(normals) <- "double"
    storage.mode(it) <- "integer"
    out <- .Fortran("adnormals",v,ncol(v),it,ncol(it),nrow(it),normals=normals,angle)$normals
    normals <- out
    normals[4,] <- 1
    x$normals <- normals
    return(x)
}
#' @export facenormals
#' @rdname adnormals
facenormals <- function(x) 
{
    barymesh <- list()
    barymesh$vb <- rbind(t(barycenter(x)),1)
    
    v <- x$vb
    ## Make sure v is homogeneous with unit w
    if (nrow(v) == 3)
        v <- rbind(v, 1)
    else
        v <- t( t(v)/v[4,] )
    v <- v[1:3,]
    it <- x$it
    normals <- it*0
    storage.mode(v) <- "double"
    storage.mode(normals) <- "double"
    storage.mode(it) <- "integer"
    out <- .Fortran("facenormals",v,ncol(v),it,ncol(it),nrow(it),normals=normals)$normals
    normals <- out
    class(barymesh) <- "mesh3d"
    barymesh$normals <- normals
    
    return(barymesh)
}


