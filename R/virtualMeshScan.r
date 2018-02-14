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
#' @examples
#' SCP1 <- file2mesh(system.file("extdata","SCP1.ply",package="Morpho"))
#' viewpoints <- read.fcsv(system.file("extdata","SCP1_Endo.fcsv",package="Morpho"))
#' visivert <- getVisibleVertices(SCP1,viewpoints)
#' @importFrom Rvcg vcgRaySearch
#' @importFrom parallel mclapply clusterExport makeCluster
#' @export
getVisibleVertices <- function(mesh,viewpoints, offset = 0.001,cores=1) {
    mesh <- mesh0 <- vcgUpdateNormals(mesh)
 ##    mesh0 <- meshOffset(mesh, offset)
    if (is.vector(viewpoints))
        if (length(viewpoints)== 3)
            viewpoints <- t(viewpoints)
    i <- 0
    parfun <- function(i) {        
        normals <- c(viewpoints[i,],0) - mesh0$vb
        mesh0$normals <- normals
        mesh0 <- meshOffset(mesh0,offset)
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
    
    if (.Platform$OS.type == "windows" && cores > 1) {
        cl <- makeCluster(cores, type='PSOCK')
        registerDoParallel(cl)
        clusterExport(cl=cl, list("mesh0", "offset", "viewpoints","meshOffset","vcgRaySearch"),envir=environment())
        outvec <- foreach(i = 1:nrow(viewpoints)) %dopar% parfun(i)
        stopCluster(cl)
        registerDoSEQ()
    } else
        outvec <- parallel::mclapply(1:nrow(viewpoints),parfun,mc.cores=cores)

    visible <- unique(unlist(outvec))
    ## invisible <- (1:Rvcg::nverts(mesh))[-visible]
    return(visible)
}

#' remove all parts of a triangular mesh, not visible from a set of viewpoints
#'
#' remove all parts of a triangular mesh, not visible from a set of viewpoints
#' @param x triangular mesh of class 'mesh3d'
#' @param viewpoints vector or k x 3  matrix containing a set of viewpoints
#' @param offset value to generate an offset at the meshes surface (see notes)
#' @param cores integer: number of cores to use (not working on windows)
#' 
#' @note The function tries to filter out all vertices where the line connecting each vertex with the viewpoints intersects with the mesh itself. As, technically speaking this always occurs at a distance of value=0, a mesh with a tiny offset is generated to avoid these false hits.
#' @return returns a list containing subsets of the original mesh
#' \item{visible}{the parts visible from at least one of the viewpoints}
#' \item{invisible}{the parts not visible from the viewpoints}
#' @examples
#' SCP1 <- file2mesh(system.file("extdata","SCP1.ply",package="Morpho"))
#' viewpoints <- read.fcsv(system.file("extdata","SCP1_Endo.fcsv",package="Morpho"))
#' ## Create a quick endocast
#' quickEndo <- virtualMeshScan(SCP1,viewpoints)
#' \dontrun{
#' rgl::shade3d(quickEndo$visible,col="orange")
#' rgl::shade3d(SCP1,col="white",alpha=0.5)
#' }
#' @importFrom Rvcg vcgRaySearch
#' @importFrom parallel mclapply
#' @export
virtualMeshScan <- function(x,viewpoints,offset=0.001,cores=1) {
    getvisible <- getVisibleVertices(mesh=x,viewpoints=viewpoints,offset=offset,cores=cores)
    visible <- rmVertex(x,getvisible,keep=TRUE)
    invisible <- rmVertex(x,getvisible,keep=FALSE)
    return(list(visible=visible,invisible=invisible))
}

#' Get viewpoints on a sphere around a 3D mesh
#'
#' Get viewpoints on a sphere around a 3D mesh to be used with virtualMeshScan
#' @param x  triangular mesh of class 'mesh3d'
#' @param n number of viewpoint to generate
#' @param inflate factor for the size of the sphere: \code{inflate=1} means that the sphere around the object just touches the point farthest away from the mesh's centroid.
#' @param radius defines a fix radius for the sphere (overrides arg \code{inflate}).
#' @param subdivision parameter passed to \code{\link{vcgSphere}}
#' @param PCA logical: if TRUE, the sphere will be deformed to match the principle axes of the mesh. NOTE: this may result in the sphere not necessarily completely enclosing the mesh.
#' @return a list containing
#' \item{viewpoints}{n x 3 matrix containing viewpoints.}
#' \item{sphere}{sphere from which the points are sampled}
#' \item{radius}{radius of the sphere}
#' @examples
#' data(boneData)
#' vp <- getOuterViewpoints(skull_0144_ch_fe.mesh,n=100)
#' \dontrun{
#' require(rgl)
#' shade3d(skull_0144_ch_fe.mesh,col="white")
#' spheres3d(vp$viewpoints)
#' wire3d(vp$sphere)
#' }
#' ### Fit to principal axes
#' vppca <- getOuterViewpoints(skull_0144_ch_fe.mesh,n=100,PCA=TRUE,inflate=1.5)
#' \dontrun{
#' require(rgl)
#' shade3d(skull_0144_ch_fe.mesh,col="white")
#' spheres3d(vppca$viewpoints)
#' wire3d(vppca$sphere)
#' }
#' @importFrom Rvcg vcgSphere vcgSample nverts
#' @export
getOuterViewpoints <- function(x,n,inflate=1.5,radius=NULL,subdivision=3,PCA=FALSE) {
    verts <- vert2points(x)
     
    mcenter <- colMeans(verts)
    mydists <- max(sqrt(rowSums((sweep(verts,2,mcenter))^2)))
    mysphere <- vcgSphere(subdivision)
    if (PCA) {
        mypca <- prcompfast(verts,retx=FALSE)
        ev <- mypca$sdev/(mypca$sdev[1])
        pcatrafo <- diag(c(ev,1))
        pcatraforot <- mypca$rotation
        pcatraforot <- rbind(cbind(pcatraforot,c(0,0,0)),c(0,0,0,1))
        mysphere <- applyTransform(mysphere,pcatraforot%*%pcatrafo%*%t(pcatraforot))
    }
    #mysphere$vb[1:3,] <- mysphere$vb[1:3,]*distance*mydists
    trafo <- diag(4)
    trafo[1:3,4] <- mcenter
    if (is.null(radius))
        radius <- mydists*inflate
    mysphere <- scalemesh(mysphere,radius,center = "none")

    mysphere <- applyTransform(mysphere,trafo)
   
        
        
    if (n > nverts(mysphere))
        stop("n exceeds the number of vertices of the sphere, increase it by setting subdivision to a higher value")
    viewpoints <- vcgSample(mysphere,n)
    if (!PCA)
        message(paste("Radius of the sphere =",radius))
    return(list(viewpoints=viewpoints,sphere=mysphere,radius=radius))
    
}
