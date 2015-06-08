#' merge multiple triangular meshes into a single one
#' 
#' merge multiple triangular meshes into a single one, preserving color and
#' vertex normals.
#' 
#' 
#' @param \dots triangular meshes of class \code{'mesh3d'} to merge or a list
#' of triangular meshes.
#' 
#' @return returns the meshes merged into a single one.
#' @seealso \code{\link{mesh2ply}, \link{file2mesh}, \link{ply2mesh}}
#' 
#' @examples
#' 
#' require(rgl)
#' data(boneData)
#' data(nose)
#' mergedMesh <- mergeMeshes(shortnose.mesh, skull_0144_ch_fe.mesh)
#' \dontrun{
#' require(rgl)
#' shade3d(mergedMesh, col=3)
#' }
#' 
#' @export
mergeMeshes <- function(...)
{
    args <- list(...)
    if (length(args) == 1 && !inherits(args[[1]],"mesh3d"))
        args <- unlist(args, recursive = FALSE)

    argc <- length(args)
        if (argc < 2)
            stop("at least two arguments needed")
    outmesh <- args[[1]]
    if (dim(outmesh$vb)[1] == 3)
        outmesh$vb <- rbind(outmesh$vb, 1)
    if (!is.null(outmesh$normals))
        if (dim(outmesh$normals)[1] == 3)
            outmesh$normals <- rbind(outmesh$normals, 1)
    if (!inherits(outmesh, "mesh3d"))
        stop("input must be meshes of class 'mesh3d'")
    
    for (i in 2 : argc) {
        if (!inherits(args[[i]], "mesh3d"))
            stop("input must be meshes of class 'mesh3d'")
        tmpmesh <- args[[i]]
        nvbout <- ncol(outmesh$vb)
        nvbtmp <- ncol(tmpmesh$vb)
        # bind vertices
        if (dim(tmpmesh$vb)[1] == 3)
            tmpmesh$vb <- rbind(tmpmesh$vb,1)
        outmesh$vb <- cbind(outmesh$vb, tmpmesh$vb)
        # bind normals
        if (!is.null(tmpmesh$normals)) {
            if (dim(tmpmesh$normals)[1] == 3)
                tmpmesh$normals <- rbind(tmpmesh$normals, 1)
        }
        if (is.null(outmesh$normals) && !is.null(tmpmesh$normals)) {
            outmesh$normals <- cbind(rbind(matrix(0, 3, nvbout), 1), tmpmesh$normals)
        } else if (!is.null(outmesh$normals) && is.null(tmpmesh$normals)) {
            outmesh$normals <- cbind(outmesh$normals, rbind(matrix(0, 3, nvbtmp), 1), tmpmesh$normals)
        } else
            outmesh$normals <- cbind(outmesh$normals, tmpmesh$normals)

        ## handle faces
        nitout <- ncol(outmesh$it)
        nittmp <- ncol(tmpmesh$it)
        if (!is.null(tmpmesh$it))
            tmpmesh$it <- tmpmesh$it+nvbout
        outmesh$it <- cbind(outmesh$it,tmpmesh$it)

        ## handle colors
        if (!is.null(tmpmesh$material$color) && is.null(outmesh$material$color) && !is.null(outmesh$it)) {
            outmesh$material$color <- matrix("#FFFFFF", 3, nitout)
        } else if (is.null(tmpmesh$material$color) && !is.null(tmpmesh$it) && !is.null(outmesh$material$color))
            tmpmesh$material$color <- matrix("#FFFFFF", 3, nittmp)
        outmesh$material$color <- cbind(outmesh$material$color, tmpmesh$material$color)
    }
    return(outmesh)
}
        
                                                            
                                     
            
