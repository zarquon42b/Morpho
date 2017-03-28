meshOffset <- function (mesh, offset) {
    if (is.null(mesh$normals)) 
        mesh <- vcgUpdateNormals(mesh)
    mesh$vb[1:3, ] <- mesh$vb[1:3, ] + t(offset * t(mesh$normals[1:3, ]))
    invisible(mesh)
}

##

#' find vertices visible from a given viewpoints
#'
#' find vertices visible from a given viewpoints
#' @param mesh triangular mesh of class 'mesh3d'
#' @param viewpoints vector or k x 3  matrix containing a set of viewpoints
#' @param offset value to generate an offset at the meshes surface (see notes)
#' @param cores integer: number of cores to use (not working on windows)
#' @note The function tries to filter out all vertices where the line connecting each vertex with the viewpoints intersects with the mesh itself. As, technically speaking this always occurs at a distance of value=0, a mesh with a tiny offset is generated to avoid these false hits.
#' @return a vector with (1-based) indices of points visible from at least one of the viewpoints
#' @importFrom Rvcg vcgRaySearch
#' @importFrom parallel mclapply
#' @export
getVisibleVertices <- function(mesh,viewpoints, offset = 0.001,cores=1) {
    mesh <- vcgUpdateNormals(mesh)
    mesh0 <- meshOffset(mesh, offset)
    if (is.vector(viewpoints))
        if (length(viewpoints)== 3)
            viewpoints <- t(viewpoints)
    parfun <- function(i) {        
        normals <- c(viewpoints[i,],0) - mesh0$vb
        mesh0$normals <- normals
        rs <- vcgRaySearch(mesh0, mesh)
        tmp <- as.logical(rs$quality)
        hitfaces <- which(tmp)
        postnormals <- c(viewpoints[i,],0) - rs$vb
        good <- angM(t(normals[,hitfaces]),t(postnormals[,hitfaces]))
        good <- which(good > pi/2)
        out <- tmp
        out[tmp] <- FALSE
        out[!tmp] <- TRUE
        out[hitfaces[good]] <- TRUE
        return(which(out))
    }
    outvec <- parallel::mclapply(1:nrow(viewpoints),parfun,mc.cores=cores)
    outvec <- unique(unlist(outvec))
    return(outvec)
}

#' remove all parts of a triangular mesh, not visible from a set of viewpoints
#'
#' remove all parts of a triangular mesh, not visible from a set of viewpoints
#' @param x triangular mesh of class 'mesh3d'
#' @param viewpoints vector or k x 3  matrix containing a set of viewpoints
#' @param offset value to generate an offset at the meshes surface (see notes)
#' @param cores integer: number of cores to use (not working on windows)
#' @note The function tries to filter out all vertices where the line connecting each vertex with the viewpoints intersects with the mesh itself. As, technically speaking this always occurs at a distance of value=0, a mesh with a tiny offset is generated to avoid these false hits.
#' @return a subset of the original mesh
#' @examples
#' SCP1 <- file2mesh(system.file("extdata","SCP1.ply",package="Morpho"))
#' viewpoints <- read.fcsv(system.file("extdata","SCP1_Endo.fcsv",package="Morpho"))
#' ## Create a quick endocast
#' quickEndo <- scanMeshFromViewpoints(SCP1,viewpoints)
#' \dontrun{
#' rgl::shade3d(quickEndo,col="orange")
#' rgl::shade3d(SCP1,col="white",alpha=0.5)
#' }
#' @importFrom Rvcg vcgRaySearch
#' @importFrom parallel mclapply
#' @export
scanMeshFromViewpoints <- function(x,viewpoints,offset=0.001,cores=1) {
    visible <- getVisibleVertices(mesh=x,viewpoints=viewpoints,offset=offset,cores=cores)
    out <- rmVertex(x,visible,keep = T)
    return(out)
}
