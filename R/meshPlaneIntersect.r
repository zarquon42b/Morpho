#' get only those faces intersecting a plane
#'
#' get only those faces intersecting a plane
#' 
#' @param v1 numeric vector of length=3 specifying a point on the separating plane
#' @param v2 numeric vector of length=3 specifying a point on the separating plane
#' @param v3 numeric vector of length=3 specifying a point on the separating plane
#' @return returns a mesh with only those faces intersecting a plane
#' @export
meshPlaneIntersect <- function(mesh, v1, v2, v3) {
    
    pointcloud <- vert2points(mesh)
    updown <- cutSpace(pointcloud, v1, v2, v3)
    ## get infos about up and down
    upface <- getFaces(mesh,which(updown))
    downface <- getFaces(mesh,which(!updown))
    nit <- 1:ncol(mesh$it)
    facesInter <- which(as.logical((nit %in% upface) * nit %in% downface))
    mesh$it <- mesh$it[,facesInter]
    mesh <- rmUnrefVertex(mesh)
    return(mesh)
}

#' find indices of faces that contain specified vertices
#'
#' find indices of faces that contain specified vertices
#' 
#' @param mesh triangular mesh of class "mesh3d"
#' @param index vector containing indices of vertices
#' @return vector of face indices
#' @export
getFaces <- function(mesh,index) {
    it <- mesh$it
    itdim <- dim(it)
    lRm <- length(index)
    vbn <- dim(mesh$vb)[2]
    indOrig <- 1:vbn
    indOut <- indOrig*0
    indNew <- 1:(vbn-lRm)     
    indOut[-index] <- indNew
    
    facefun <- function(x)
        {
            x <- indOut[x]
            return(x)
        }
    if (!is.null(it)) {
        it <- matrix(facefun(it),itdim)
        checkface <- .Call("face_zero",it)
                                        #checkface <- .Fortran("face_zero",it,itdim[2],checkface)[[3]]
        invalface <- which(checkface == 0)
        return(invalface)
    } else {
        stop("mesh contains no triangular faces")
    }
        
}
