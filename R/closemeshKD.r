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
#' \item{quality }{vector: containing distances to target}
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
#' @keywords ~kwd1 ~kwd2
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
            mesh <- updateNormals(mesh)
        if (is.matrix(x)) {
            matr <- x
            x <- list()
            x$vb <- rbind(t(matr),1)
        } else {
            matr <- t(x$vb[1:3,])
        }
        vb <- (mesh$vb[1:3,])
        it <- (mesh$it)
        nvb <- ncol(vb)
        nit <- ncol(it)        
        nmat <- nrow(matr)
        dist <- rep(0,nmat)
        fptr <- dist
        bary <- barycenter(mesh)
        clostInd <- mcNNindex(bary,matr,k=k,cores=cores,...)
        storage.mode(k) <- "integer"
        clost <- c(0,0,0)
        storage.mode(fptr) <- "integer"
        storage.mode(it) <- "integer"
        storage.mode(nvb) <- "integer"    
        storage.mode(nit) <- "integer"
        storage.mode(method) <- "integer"
        storage.mode(nmat) <- "integer"
        storage.mode(matr) <- "double"
        storage.mode(vb) <- "double"
        storage.mode(dist) <- "double"
        storage.mode(clost) <- "double"
        storage.mode(clostInd) <- "integer"
        if (barycoords) {
            baryc <- t(matr)
            baryn <- nmat
        } else {
            baryc <- c(0,0,0)
            baryn <- as.integer(1)
        }
        outmatr <- matr
        region <- fptr
        out <- .Fortran("matr_meshKD",ioMat=matr,nmat,vb,nvb,it,nit,clostInd,k,dist=dist,faceptr=fptr,region,normals=mesh$normals[1:3,],sign,outnormals=t(matr),method, barycoords=baryc, baryn)
        gc()
        x$vb[1:3,] <- t(out$ioMat)
        x$quality <- out$dist
        x$normals <- rbind(out$outnormals,1)
        x$faceptr <- out$faceptr
        if (barycoords)
            x$barycoords <- out$barycoords
        class(x) <- "mesh3d"
        return(x)
    }
