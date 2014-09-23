#' Project coordinates onto a target triangular surface mesh.
#' 
#' For a set of 3D-coordinates the closest matches on a target surface are
#' determined and normals at as well as distances to that point are calculated.
#' 
#' The search for the clostest point is designed as follows: Calculate the
#' barycenter of each target face. For each coordinate of x, determine the k
#' closest barycenters and calculate the distances to the closest point on
#' these faces.
#' 
#' @param x k x 3 matrix containing 3D-coordinates or object of class
#' \code{mesh3d}.
#' 
#' @param mesh triangular surface mesh stored as object of class \code{mesh3d}.
#' @param k neighbourhood of kd-tree to search - the larger, the slower - but
#' the more likely the absolutely closest point is hit.
#' @param sign logical: if TRUE, signed distances are returned.
#' @param barycoords logical: if \code{TRUE}, barycentric coordinates of the
#' hit points are returned.
#' @param cores integer: how many cores to use for the search algorithm.
#' @param method integer: either 0 or 1, if 0 ordinary Euclidean distance is
#' used, if 1, the distance suggested by Moshfeghi(1994) is calculated.
#' @param \dots additional arguments. currently unavailable.
#' @return returns an object of class \code{mesh3d}.  with:
#' \item{vb }{4xn matrix containing n vertices as homolougous coordinates}
#' \item{normals }{4xn matrix containing vertex normals}
#' \item{quality }{vector: containing distances to target. In case of \code{method=1}, this is not the Euclidean distance but the distance of the reference point to the faceplane (orthogonally projected) plus the distance to the closest point on one of the face's edges (the target point). See the literature cited below for details.}
#' \item{it }{4xm matrix containing vertex indices forming triangular faces.Only available, when x is a mesh}
#' @author Stefan Schlager
#' @seealso \code{\link{ply2mesh}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling.
#' 
#' Moshfeghi M, Ranganath S, Nawyn K. 1994. Three-dimensional elastic matching
#' of volumes IEEE Transactions on Image Processing: A Publication of the IEEE
#' Signal Processing Society 3:128-138.
#' 
#' @examples
#' 
#' require(rgl)
#' data(nose)
#' out <- closemeshKD(longnose.lm,shortnose.mesh,sign=TRUE)
#' ### show distances - they are very small because
#' ###longnose.lm is scaled to unit centroid size.
#' hist(out$quality)
#' 
#' @export
closemeshKD <- function(x, mesh, k=50, sign=FALSE, barycoords=FALSE, cores=1, method=0,...)
    {
        if(.Platform$OS.type == "windows")
            cores <- 1
        if (is.null(mesh$normals))
            mesh <- vcgUpdateNormals(mesh)
        if (is.matrix(x)) {
            matr <- x
            x <- list()
            x$vb <- rbind(t(matr),1)
        } else {
            matr <- t(x$vb[1:3,])
        }
        vb <- (mesh$vb[1:3,])
        if (!is.null(mesh$it))
            it <- mesh$it-1
        else
            stop("mesh has no triangular faces")
        
        bary <- barycenter(mesh)
        clostInd <- mcNNindex(bary,matr,k=k,cores=cores,...)-1
        normals <- mesh$normals[1:3,]
        out <- .Call("points2mesh", t(matr), vb, it, normals, t(clostInd), sign, barycoords, method)
        gc()
        x$vb[1:3,] <- out$clost
        x$quality <- out$dists
        x$normals <- rbind(out$normals,1)
        x$faceptr <- out$faceptr+1
        if (barycoords)
            x$barycoords <- out$barycoords
        class(x) <- "mesh3d"
        return(x)
    }
