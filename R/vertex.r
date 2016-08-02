#' some little helpers for vertex operations on triangular meshes
#' 
#' extract vertex coordinates from meshes, find and/or remove (unreferenced)
#' vertices from triangular meshes
#' 
#' \code{unrefVertex} finds unreferenced vertices in triangular meshes of class
#' \code{mesh3d}.
#' 
#' \code{rmVertex} removes specified vertices from triangular meshes.
#' 
#' \code{vert2points} extacts vertex coordinates from triangular meshes.
#' 
#' \code{rmUnrefVertex} removes unreferenced vertices from triangular meshes.
#' @title  some little helpers for vertex operations on triangular meshes
#' @param mesh triangular mesh of class \code{mesh3d}.
#' @param index vector containing indices of vertices to be removed.
#' @param keep logical: if \code{TRUE}, the vertices specified by \code{index}
#' are kept and the rest is removed.
#' @param silent logical: suppress output about info on removed vertices.
#' @return \code{unrefVertex}: vector with indices of unreferenced vertices.
#' 
#' \code{rmVertex}: returns mesh with specified vertices removed and faces and
#' normals updated.
#' 
#' \code{vert2points}: k x 3 matrix containing vertex coordinates.
#' 
#' \code{rmUnrefVertex}: mesh with unreferenced vertices removed.
#' @author Stefan Schlager
#' @seealso \code{\link{ply2mesh}}, \code{\link{file2mesh}}
#' 
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' testmesh <- rmVertex(shortnose.mesh,1:50) ## remove first 50 vertices
#' \dontrun{
#' shade3d(testmesh,col=3)### view result
#' }
#' testmesh$vb <- cbind(testmesh$vb,shortnose.mesh$vb[,1:50]) ## add some unreferenced vertices
#' \dontrun{
#' points3d(vert2points(testmesh),col=2)## see the vertices in the holes?
#' }
#' cleanmesh <- rmUnrefVertex(testmesh)## remove those lonely vertices!
#' \dontrun{
#' rgl.pop()
#' points3d(vert2points(cleanmesh),col=2) ### now the holes are empty!!
#' }
#' 
#' @rdname vertex
#' @export
unrefVertex <- function(mesh)
    {
        it <- mesh$it
        vind <- 1:dim(mesh$vb)[2]
        unref <- which(! vind %in% it)
        return(unref)
    }
#' @rdname vertex
#' @export
rmVertex <- function(mesh,index,keep=FALSE) {
    if (! keep) {
        index <- unique(index)
        it <- mesh$it
        itdim <- dim(it)
        lRm <- length(index)
        vbn <- dim(mesh$vb)[2]
        indOrig <- 1:vbn
        indOut <- indOrig*0
        indNew <- 1:(vbn-lRm)     
        indOut[-index] <- indNew
        
        facefun <- function(x) {
            x <- indOut[x]
            return(x)
        }
        if (!is.null(it)) {
            it <- matrix(facefun(it),itdim)
            checkface <- .Call("face_zero",it)
                                        #checkface <- .Fortran("face_zero",it,itdim[2],checkface)[[3]]
            invalface <- which(checkface == 0) 
            if (length(invalface) > 0) {
                if (length(invalface) == ncol(it)) {
                    mesh$material <- NULL
                    mesh$it <- NULL
                } else
                    mesh$it <- it[,-invalface]
                if(!is.null(mesh$material$color))
                    mesh$material$color <- mesh$material$color[,-invalface]
            } else {
                mesh$it <- it
            }
            
            if (0 %in% dim(it))
                mesh$it <- NULL
        }
       
        mesh$vb <- mesh$vb[,-index]
        if (!is.null(mesh$it))
            mesh <- vcgUpdateNormals(mesh)
        else
            mesh$normals <- NULL
    } else {
        mesh <- rmVertex(mesh,c(1:ncol(mesh$vb))[-unique(index)],keep = F)
    }
    return(mesh)

}
#' @rdname vertex
#' @export
vert2points <- function(mesh)
    {
        out <- t(mesh$vb[1:3,])
        return(out)
    }
#' @rdname vertex
#' @export
rmUnrefVertex <- function(mesh, silent=FALSE)
    {
        unref <- unrefVertex(mesh)
        lunr <- length(unref)
        if (lunr > 0) {
            mesh <- rmVertex(mesh,unref)
            if (!silent)
                cat(paste(" removed",lunr,"unreferenced vertices\n"))
        } else {
            if (!silent)
                cat(" nothing to remove\n")
        }
        return(mesh)
    }

