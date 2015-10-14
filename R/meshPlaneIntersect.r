#' get intersections between mesh and a plane
#'
#' get intersections between mesh and a plane
#' @param mesh triangular mesh of class "mesh3d"
#' @param v1 numeric vector of length=3 specifying a point on the separating plane
#' @param v2 numeric vector of length=3 specifying a point on the separating plane
#' @param v3 numeric vector of length=3 specifying a point on the separating plane
#' @param normal plane normal (overrides specification by v2 and v3)
#' @return returns the intersections of edges and the plane
#' @examples
#' data(nose)
#' v1 <- shortnose.lm[1,]
#' v2 <- shortnose.lm[2,]
#' v3 <- shortnose.lm[3,]
#' intersect <- meshPlaneIntersect(shortnose.mesh,v1,v2,v3)
#' \dontrun{
#' require(rgl)
#' wire3d(shortnose.mesh)
#' spheres3d(shortnose.lm[1:3,],col=2)#the plane
#' spheres3d(intersect,col=3,radius = 0.2)#intersections
#' }
#' @importFrom Rvcg vcgGetEdge
#' @export
meshPlaneIntersect <- function(mesh, v1, v2=NULL, v3=NULL,normal=NULL) {
    
    pointcloud <- vert2points(mesh)
    updown <- cutSpace(pointcloud,v1=v1,v2=v2,v3=v3,normal=normal)
    ## get infos about up and down
    upface <- getFaces(mesh,which(updown))
    downface <- getFaces(mesh,which(!updown))
    nit <- 1:ncol(mesh$it)
    facesInter <- which(as.logical((nit %in% upface) * nit %in% downface))
    mesh$it <- mesh$it[,facesInter]
    mesh <- rmUnrefVertex(mesh,silent=TRUE)
    edges <- as.matrix(vcgGetEdge(mesh)[,1:2])
    pointcloud <- vert2points(mesh)
    out <- edgePlaneIntersect(pointcloud,edges,v1=v1,v2=v2,v3=v3,normal=normal)
    return(out)
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
    index <- unique(index)
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
        return(numeric(0))
    }
        
}

edgePlaneIntersect <- function(pointcloud,edges,v1, v2=NULL, v3=NULL,normal=NULL) {
    if (!is.null(normal)) {
        tangent <- tangentPlane(normal)
        v2 <- v1+tangent$z
        v3 <- v1+tangent$y
    }
    e1 <- v2-v1
    e2 <- v3-v1
    e1 <- e1/sqrt(sum(e1^2))
    e2 <- e2/sqrt(sum(e2^2))
    normal <- crossProduct(e1,e2)
    normal <- normal/sqrt(sum(normal^2))
    e2a <- crossProduct(e1,normal)
    e2a <- e2a/sqrt(sum(e2a^2))
    Ep <- cbind(e1,e2a)
    edges <- edges-1
    pointcloud0 <- sweep(pointcloud,2,v1)
    orthopro <- t(Ep%*%t(Ep)%*%t(pointcloud0))    
    diff <- orthopro-pointcloud0
    out <- .Call("edgePlane",pointcloud,diff,edges)
    return(out)
}

binfun <- function(coords,bs=1,margin=10){
    pca <- prcomp(coords)
    S0 <- pca$x[,1:2]
    i <- floor((S0[,2])/bs)
    j <- floor(S0[,1]/bs)
    irange <- range(i)
    jrange <- range(j)
    i <- i-irange[1]+1+margin
    j <- j-jrange[1]+1+margin
    irange <- range(i)
    jrange <- range(j)
    image <- matrix(0,irange[2]+margin,jrange[2]+margin)
    ij <- unique(cbind(i,j))
    for (i in 1:nrow(ij))
        image[ij[i,1],ij[i,2]] <- 1
    return(image)
}
